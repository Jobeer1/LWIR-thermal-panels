"""
geometry.py — 3-D cavity geometry definitions.

Three geometry types are provided:

  RectPit3D       — Rectangular trench, width × depth × height (µm).
                    Models a slit trench (kept for testing/legacy).

  FrustumCavity3D — Tapered cylindrical frustum (truncated cone) with
                    separate base and top radii.

  CNTForestCell   — Unit cell of a vertically-aligned CNT forest.
                    Square pitch, exterior of solid frustum pillar.
                    Periodic boundary conditions at cell walls.
                    Correct 3-D areas for cavity enhancement factor.

  HoneycombCavityCell — Unit cell of a hexagonally-packed cylindrical
                    cavity panel. The cavity is a cylinder bored INTO a
                    substrate. Walls and base are the high-emissivity
                    coating. Aperture is the circular opening at top.
                    Uses hex close-packing to scale to a panel.

All geometries expose:
  .area_walls     : total emitting wall area (m² per unit cell)
  .area_base      : base area (m² per unit cell)
  .area_aperture  : aperture area (m² per unit cell)
  .height         : cavity depth (m)
  .kind           : string identifier

Ray intersection is provided by geometry-specific methods used by ray_tracer.py.
"""

import math
import numpy as np
from dataclasses import dataclass
from typing import Optional, Tuple


# ---------------------------------------------------------------------------
# Helper: ray-axis-aligned-box intersection (used for rect pit walls)
# ---------------------------------------------------------------------------

def _aabb_intersect(pos: np.ndarray, direction: np.ndarray,
                    lo: np.ndarray, hi: np.ndarray
                    ) -> Tuple[float, Optional[np.ndarray]]:
    """Ray vs axis-aligned bounding box.  Returns (t, normal) for closest
    positive-t face hit, or (inf, None) if no hit.  Skips the top face
    (hi[2]) — that is the aperture handled separately by the tracer.
    """
    EPS = 1e-12
    t_min = 0.0
    t_max = float('inf')
    hit_normal = None

    for i in range(3):
        if abs(direction[i]) < EPS:
            if pos[i] < lo[i] or pos[i] > hi[i]:
                return float('inf'), None
        else:
            t1 = (lo[i] - pos[i]) / direction[i]
            t2 = (hi[i] - pos[i]) / direction[i]
            sign = 1 if direction[i] > 0 else -1
            n1 = np.zeros(3); n1[i] = -sign
            n2 = np.zeros(3); n2[i] = sign
            if t1 > t2:
                t1, t2 = t2, t1
                n1, n2 = n2, n1
            if t1 > t_min:
                t_min = t1
                hit_normal = n1
            t_max = min(t_max, t2)
            if t_min > t_max:
                return float('inf'), None

    if t_min <= 0.0:
        return float('inf'), None
    return t_min, hit_normal


# ---------------------------------------------------------------------------
# RectPit3D
# ---------------------------------------------------------------------------

@dataclass
class RectPit3D:
    """Rectangular trench (µm inputs, internal storage in metres).

    The aperture is the full width × depth face at z = height.
    Walls are the four vertical sides.

    Parameters
    ----------
    width_um  : trench width in µm  (x-axis)
    depth_um  : trench depth into page in µm  (y-axis)
    height_um : trench height in µm  (z-axis, aperture at top)
    """
    width_um: float
    depth_um: float
    height_um: float
    kind: str = 'rect_pit'

    def __post_init__(self):
        self.W = self.width_um * 1e-6
        self.D = self.depth_um * 1e-6
        self.H = self.height_um * 1e-6
        self.area_aperture = self.W * self.D
        self.area_base     = self.W * self.D
        self.area_walls    = 2.0 * (self.W + self.D) * self.H
        self.height        = self.H

        # Waveguide modal cutoff (peer-review physics): the lowest rectangular
        # mode is limited by the narrower transverse dimension.
        self.lambda_c_um = 2.0 * min(self.width_um, self.depth_um)
        self.lambda_c    = self.lambda_c_um * 1e-6

    def channel_cutoff_wavelength_um(self) -> float:
        """Modal cutoff wavelength (µm). Emission λ > λ_c is evanescent."""
        return self.lambda_c_um

    @property
    def cavity_enhancement(self) -> float:
        """Geometric ratio of (wall + base) area to aperture area."""
        return (self.area_walls + self.area_base) / self.area_aperture

    def sample_point_on_walls(self) -> Tuple[np.ndarray, np.ndarray]:
        """Random point on one of the four side walls; returns (pos_m, outward_normal)."""
        # perimeter weights
        px = self.W / (2 * (self.W + self.D))
        py = self.D / (2 * (self.W + self.D))
        r = np.random.random()
        z = np.random.uniform(0.0, self.H)
        if r < px:         # x=0 wall
            return np.array([0.0, np.random.uniform(0.0, self.D), z]), np.array([1.0, 0.0, 0.0])
        elif r < 2*px:     # x=W wall
            return np.array([self.W, np.random.uniform(0.0, self.D), z]), np.array([-1.0, 0.0, 0.0])
        elif r < 2*px+py:  # y=0 wall
            return np.array([np.random.uniform(0.0, self.W), 0.0, z]), np.array([0.0, 1.0, 0.0])
        else:              # y=D wall
            return np.array([np.random.uniform(0.0, self.W), self.D, z]), np.array([0.0, -1.0, 0.0])

    def sample_point_on_base(self) -> Tuple[np.ndarray, np.ndarray]:
        """Random point on base (z=0); returns (pos_m, outward_normal)."""
        return (np.array([np.random.uniform(0.0, self.W),
                          np.random.uniform(0.0, self.D), 0.0]),
                np.array([0.0, 0.0, 1.0]))

    def next_hit(self, pos: np.ndarray, direction: np.ndarray
                 ) -> Tuple[float, str, np.ndarray]:
        """Find next ray-boundary intersection inside the pit.

        Returns (t, surface_type, outward_normal).
        surface_type ∈ {'aperture', 'wall', 'base'}
        """
        EPS = 1e-14
        # Aperture (top): z = H
        t_top = (self.H - pos[2]) / direction[2] if direction[2] > EPS else float('inf')
        # Base (bottom): z = 0
        t_bot = (0.0 - pos[2]) / direction[2] if direction[2] < -EPS else float('inf')
        # Side walls
        t_xm = (0.0 - pos[0]) / direction[0] if direction[0] < -EPS else float('inf')
        t_xp = (self.W - pos[0]) / direction[0] if direction[0] > EPS else float('inf')
        t_ym = (0.0 - pos[1]) / direction[1] if direction[1] < -EPS else float('inf')
        t_yp = (self.D - pos[1]) / direction[1] if direction[1] > EPS else float('inf')

        # Build (t, surface, normal) candidates — skip non-positive
        candidates = []
        if t_top > EPS:   candidates.append((t_top, 'aperture', np.array([0., 0., 1.])))
        if t_bot > EPS:   candidates.append((t_bot, 'base',     np.array([0., 0., 1.])))
        if t_xm > EPS:    candidates.append((t_xm,  'wall',     np.array([1., 0., 0.])))
        if t_xp > EPS:    candidates.append((t_xp,  'wall',     np.array([-1.,0., 0.])))
        if t_ym > EPS:    candidates.append((t_ym,  'wall',     np.array([0., 1., 0.])))
        if t_yp > EPS:    candidates.append((t_yp,  'wall',     np.array([0.,-1., 0.])))

        if not candidates:
            return float('inf'), 'none', np.zeros(3)
        return min(candidates, key=lambda c: c[0])


# ---------------------------------------------------------------------------
# FrustumCavity3D  (tapered / conical CNT geometry)
# ---------------------------------------------------------------------------

@dataclass
class FrustumCavity3D:
    """Tapered cylindrical frustum (truncated cone) cavity.

    The cone axis is aligned with z.  The aperture is the top face at z=H
    with radius r_top.  The base is at z=0 with radius r_base.

    Parameters (all in µm, stored internally in metres)
    ----------
    r_base_um : radius at z=0 (base, larger end for typical CNTs)
    r_top_um  : radius at z=H (aperture, smaller / tip end)
    height_um : depth of the frustum
    """
    r_base_um: float
    r_top_um:  float
    height_um: float
    kind: str = 'frustum'

    def __post_init__(self):
        self.r_base = self.r_base_um * 1e-6
        self.r_top  = self.r_top_um  * 1e-6
        self.H      = self.height_um * 1e-6
        self.height = self.H
        r_avg = (self.r_base + self.r_top) / 2.0
        slant = math.hypot(self.H, self.r_base - self.r_top)
        self.area_walls    = math.pi * (self.r_base + self.r_top) * slant
        self.area_base     = math.pi * self.r_base ** 2
        self.area_aperture = math.pi * self.r_top  ** 2
        # For a frustum the effective aperture used in epsilon_B is the top circle
        self._dr_dz = (self.r_top - self.r_base) / self.H   # dr/dz (negative if tapering in)

        # Waveguide modal cutoff (Bug #13):
        #   To ESCAPE, a mode must survive up to the narrowest end (the aperture,
        #   radius r_top), so the cutoff is set by the aperture diameter:
        #   TE11 mode of a circular guide:  λ_c = 1.706 · d_top.
        self.lambda_c_um = 1.706 * (2.0 * self.r_top_um)
        self.lambda_c    = self.lambda_c_um * 1e-6

    def channel_cutoff_wavelength_um(self) -> float:
        """Modal cutoff wavelength (µm). Emission λ > λ_c is evanescent."""
        return self.lambda_c_um

    @property
    def cavity_enhancement(self) -> float:
        return (self.area_walls + self.area_base) / self.area_aperture

    def _radius_at_z(self, z: float) -> float:
        return self.r_base + self._dr_dz * z

    def next_hit(self, pos: np.ndarray, direction: np.ndarray
                 ) -> Tuple[float, str, np.ndarray]:
        """Ray vs frustum intersection.

        The frustum wall can be described as a cone:
            x² + y² = (r_base + dr_dz * z)²
        Expanding: (dx²+dy²)t² + 2(px*dx + py*dy - dr_dz²*(pz-z0)*dz - r_base*dr_dz*dz)t + ...
        Uses standard analytic quadratic.
        Returns (t, surface_type, outward_normal).
        """
        EPS = 1e-14
        # Aperture
        t_top = (self.H - pos[2]) / direction[2] if direction[2] > EPS else float('inf')
        # Base
        t_bot = (0.0 - pos[2]) / direction[2] if direction[2] < -EPS else float('inf')

        # Cone wall: x² + y² = (r_base + k*z)², k = dr_dz
        k = self._dr_dz
        px, py, pz = pos[0], pos[1], pos[2]
        dx, dy, dz = direction[0], direction[1], direction[2]
        # Substitute x = px+dx*t etc.
        A = dx*dx + dy*dy - k*k*dz*dz
        rz = self.r_base + k * pz
        B = 2.0 * (px*dx + py*dy - k*rz*dz)
        C = px*px + py*py - rz*rz
        disc = B*B - 4.0*A*C

        cone_candidates = []
        if abs(A) < 1e-20:
            # Linear case
            if abs(B) > 1e-20:
                t_c = -C / B
                if t_c > EPS:
                    cone_candidates.append(t_c)
        elif disc >= 0.0:
            sq = math.sqrt(disc)
            for sign in (-1, 1):
                t_c = (-B + sign * sq) / (2.0 * A)
                if t_c > EPS:
                    z_hit = pz + dz * t_c
                    if 0.0 <= z_hit <= self.H:
                        cone_candidates.append(t_c)

        candidates = []
        if t_top > EPS:
            # Check inside aperture circle
            xh = pos[0] + direction[0] * t_top
            yh = pos[1] + direction[1] * t_top
            if xh*xh + yh*yh <= self.r_top**2:
                candidates.append((t_top, 'aperture', np.array([0., 0., 1.])))
        if t_bot > EPS:
            xh = pos[0] + direction[0] * t_bot
            yh = pos[1] + direction[1] * t_bot
            if xh*xh + yh*yh <= self.r_base**2:
                candidates.append((t_bot, 'base', np.array([0., 0., 1.])))
        for t_c in cone_candidates:
            xh = pos[0] + direction[0] * t_c
            yh = pos[1] + direction[1] * t_c
            zh = pos[2] + direction[2] * t_c
            rh = self._radius_at_z(zh)
            # Outward normal on cone surface (pointing inward towards axis is -n_out)
            # n_out in (x,y) plane perpendicular to cone surface
            nxy = np.array([xh, yh, 0.0])
            if np.linalg.norm(nxy) > 0:
                nxy /= np.linalg.norm(nxy)
            # Cone half-angle contribution
            n_out = np.array([nxy[0], nxy[1], -k / math.hypot(1.0, k)])
            n_out /= np.linalg.norm(n_out)
            # Ensure it points outward (away from axis)
            if np.dot(n_out[:2], np.array([xh, yh])) < 0:
                n_out = -n_out
            candidates.append((t_c, 'wall', n_out))

        if not candidates:
            return float('inf'), 'none', np.zeros(3)
        return min(candidates, key=lambda c: c[0])

    def sample_point_on_walls(self) -> Tuple[np.ndarray, np.ndarray]:
        """Uniform random point on the frustum lateral surface."""
        # Use area-proportional z sampling for a linearly tapered cone
        # Area element: dA = pi*(r_base + k*z)*slant_per_unit * dz
        k = self._dr_dz
        # Sample z proportional to r(z) via rejection
        while True:
            z = np.random.uniform(0.0, self.H)
            rz = self._radius_at_z(z)
            if np.random.random() < rz / max(self.r_base, self.r_top):
                break
        phi = np.random.uniform(0.0, 2 * math.pi)
        rz = self._radius_at_z(z)
        pos = np.array([rz * math.cos(phi), rz * math.sin(phi), z])
        # Outward normal
        nxy = np.array([math.cos(phi), math.sin(phi), 0.0])
        n_out = np.array([nxy[0], nxy[1], -k / math.hypot(1.0, k)])
        n_out /= np.linalg.norm(n_out)
        if np.dot(n_out[:2], pos[:2]) < 0:
            n_out = -n_out
        return pos, n_out

    def sample_point_on_base(self) -> Tuple[np.ndarray, np.ndarray]:
        """Uniform random point on base disk."""
        r = self.r_base * math.sqrt(np.random.random())
        phi = np.random.uniform(0.0, 2 * math.pi)
        return (np.array([r * math.cos(phi), r * math.sin(phi), 0.0]),
                np.array([0.0, 0.0, 1.0]))


# ---------------------------------------------------------------------------
# CNTForestCell
# ---------------------------------------------------------------------------

@dataclass
class CNTForestCell:
    """Unit cell of a vertically-aligned CNT forest (square pitch).

    Models the interstitial space around a solid frustum (CNT) in a 
    square periodic unit cell.
    """
    pitch_um:    float
    dia_base_nm: float
    dia_top_nm:  float
    height_um:   float
    kind: str = 'cnt_forest'

    def __post_init__(self):
        self.P  = self.pitch_um * 1e-6
        self.rb = (self.dia_base_nm * 1e-9) / 2.0
        self.rt = (self.dia_top_nm  * 1e-9) / 2.0
        self.H  = self.height_um * 1e-6
        self.height = self.H

        self._dr_dz = (self.rt - self.rb) / self.H if self.H > 0 else 0.0

        # Areas
        self.slant = math.hypot(self.H, self.rb - self.rt)
        self.area_lateral = math.pi * (self.rb + self.rt) * self.slant
        self.area_top_cap = math.pi * self.rt**2
        self.area_walls   = self.area_lateral + self.area_top_cap
        
        self.area_base    = self.P**2 - math.pi * self.rb**2
        self.area_aperture = self.P**2

        # Interstitial "channel" width between adjacent tubes (µm).
        # A tapered CNT narrows from base (dia_base) to tip (dia_top); use the
        # volume-weighted mean diameter for the characteristic gap.
        mean_diameter_um = ((self.dia_base_nm + self.dia_top_nm) / 2.0) / 1000.0
        self.gap_um = max(self.pitch_um - mean_diameter_um, 0.1 * self.pitch_um)

        # Waveguide modal cutoff (Bug #13):
        #   The interstitial space ≈ square waveguide of side gap → TE10 cutoff
        #   λ_c = 2·gap.  For sub-wavelength inter-tube gaps (typical CNT forest
        #   thermal IR) essentially ALL thermal emission is evanescent and the
        #   forest interior is a strong-mode-suppressed emitter.
        self.lambda_c_um = 2.0 * self.gap_um
        self.lambda_c    = self.lambda_c_um * 1e-6

    def channel_cutoff_wavelength_um(self) -> float:
        """Modal cutoff wavelength (µm). Emission λ > λ_c is evanescent."""
        return self.lambda_c_um

    @property
    def cavity_enhancement(self) -> float:
        return (self.area_walls + self.area_base) / self.area_aperture

    def _radius_at_z(self, z: float) -> float:
        return self.rb + self._dr_dz * z

    def next_hit(self, pos: np.ndarray, direction: np.ndarray) -> Tuple[float, str, np.ndarray]:
        EPS = 1e-12
        candidates = []

        # 1. Z-bounds (Aperture / Base / Top Cap)
        if direction[2] > EPS:
            t_z = (self.H - pos[2]) / direction[2]
            if t_z > EPS:
                # Escaping through the top
                candidates.append((t_z, 'aperture', np.array([0., 0., 1.])))
        elif direction[2] < -EPS:
            t_z = (0.0 - pos[2]) / direction[2]
            if t_z > EPS:
                # Hitting the floor
                candidates.append((t_z, 'base', np.array([0., 0., 1.])))
                
            # Check if hitting the top cap from above (e.g., incident rays)
            t_topcap = (self.H - pos[2]) / direction[2]
            if t_topcap > EPS:
                x_hit = pos[0] + direction[0] * t_topcap
                y_hit = pos[1] + direction[1] * t_topcap
                if x_hit*x_hit + y_hit*y_hit <= self.rt**2:
                    candidates.append((t_topcap, 'top_cap', np.array([0., 0., 1.])))

        # 2. Periodic Boundaries (x = ±P/2, y = ±P/2)
        if direction[0] > EPS:
            t_x = (self.P/2 - pos[0]) / direction[0]
            if t_x > EPS: candidates.append((t_x, 'periodic_x', np.array([-1., 0., 0.])))
        elif direction[0] < -EPS:
            t_x = (-self.P/2 - pos[0]) / direction[0]
            if t_x > EPS: candidates.append((t_x, 'periodic_x', np.array([1., 0., 0.])))

        if direction[1] > EPS:
            t_y = (self.P/2 - pos[1]) / direction[1]
            if t_y > EPS: candidates.append((t_y, 'periodic_y', np.array([0., -1., 0.])))
        elif direction[1] < -EPS:
            t_y = (-self.P/2 - pos[1]) / direction[1]
            if t_y > EPS: candidates.append((t_y, 'periodic_y', np.array([0., 1., 0.])))

        # 3. Frustum Exterior (Cylinder wall)
        # We are OUTSIDE the cylinder looking IN, so we want the FIRST positive root.
        k = self._dr_dz
        px, py, pz = pos[0], pos[1], pos[2]
        dx, dy, dz = direction[0], direction[1], direction[2]
        
        A = dx*dx + dy*dy - k*k*dz*dz
        rz = self.rb + k * pz
        B = 2.0 * (px*dx + py*dy - k*rz*dz)
        C = px*px + py*py - rz*rz
        disc = B*B - 4.0*A*C

        if abs(A) < 1e-20:
            if abs(B) > 1e-20:
                t_c = -C / B
                if t_c > EPS:
                    z_hit = pz + dz * t_c
                    if 0.0 <= z_hit <= self.H:
                        # Normal vector points OUT of the cylinder (into our interstitial space)
                        r_hit = self._radius_at_z(z_hit)
                        x_hit = px + dx * t_c
                        y_hit = py + dy * t_c
                        nx = x_hit / r_hit if r_hit > 0 else 0
                        ny = y_hit / r_hit if r_hit > 0 else 0
                        nz = -k
                        n_out = np.array([nx, ny, nz])
                        n_out /= np.linalg.norm(n_out)
                        candidates.append((t_c, 'wall', n_out))
        elif disc >= 0.0:
            sq = math.sqrt(disc)
            for sign in (-1, 1):
                t_c = (-B + sign * sq) / (2.0 * A)
                if t_c > EPS:
                    z_hit = pz + dz * t_c
                    if 0.0 <= z_hit <= self.H:
                        r_hit = self._radius_at_z(z_hit)
                        x_hit = px + dx * t_c
                        y_hit = py + dy * t_c
                        nx = x_hit / r_hit if r_hit > 0 else 0
                        ny = y_hit / r_hit if r_hit > 0 else 0
                        nz = -k
                        n_out = np.array([nx, ny, nz])
                        n_out /= np.linalg.norm(n_out)
                        candidates.append((t_c, 'wall', n_out))

        if not candidates:
            return float('inf'), 'none', np.zeros(3)
        return min(candidates, key=lambda c: c[0])

    def sample_point_on_walls(self) -> Tuple[np.ndarray, np.ndarray]:
        """Sample on the frustum lateral exterior OR the top cap."""
        if np.random.random() < self.area_top_cap / self.area_walls:
            # Sample on top cap
            r = self.rt * math.sqrt(np.random.random())
            phi = np.random.uniform(0.0, 2 * math.pi)
            return (np.array([r * math.cos(phi), r * math.sin(phi), self.H]),
                    np.array([0.0, 0.0, 1.0]))
        
        # Sample on lateral wall
        z = np.random.uniform(0.0, self.H)
        r = self._radius_at_z(z)
        phi = np.random.uniform(0.0, 2 * math.pi)
        pos = np.array([r * math.cos(phi), r * math.sin(phi), z])
        nx = math.cos(phi)
        ny = math.sin(phi)
        nz = -self._dr_dz
        n_out = np.array([nx, ny, nz])
        n_out /= np.linalg.norm(n_out)
        return pos, n_out

    def sample_point_on_base(self) -> Tuple[np.ndarray, np.ndarray]:
        """Uniform random point on substrate (outside tube footprint)."""
        while True:
            x = np.random.uniform(-self.P / 2, self.P / 2)
            y = np.random.uniform(-self.P / 2, self.P / 2)
            if x*x + y*y > self.rb**2:
                return np.array([x, y, 0.0]), np.array([0.0, 0.0, 1.0])


# ---------------------------------------------------------------------------
# HoneycombCavityCell
# ---------------------------------------------------------------------------

@dataclass
class HoneycombCavityCell:
    """Unit cell of a hexagonally-packed cylindrical cavity panel.

    The unit cell is a circular cylindrical hole of diameter D_um bored
    vertically into a substrate.  The aperture (open end) is at z = H.
    The base is the flat bottom at z = 0.
    Walls are the cylindrical interior surface.

    For scaling to a full panel:
        packing_fraction f:  fraction of panel area covered by apertures
                             (default 0.9069 for close-packed hexagons)
        eps_panel = f * alpha_eff + (1-f) * eps_flat

    The pitch P used internally equals D / sqrt(f * pi / (2*sqrt(3)))
    but the aperture sampler uses a square cell of side P for simplicity.
    """
    diameter_um:      float   # cavity opening diameter (µm)
    height_um:        float   # cavity depth (µm)
    wall_emissivity:  float = 0.95   # stored for reference; passed to ray_tracer separately
    packing_fraction: float = 0.9069 # hex close-pack
    wall_material:    str   = 'alumina'  # Fix 1: material for complex dielectric cutoff
    kind: str = 'honeycomb'

    def __post_init__(self):
        self.D  = self.diameter_um * 1e-6
        self.R  = self.D / 2.0
        self.H  = self.height_um * 1e-6
        self.height = self.H

        # Equivalent square cell side for aperture sampling
        # area_cell = pi * R^2 / packing_fraction
        import math as _math
        self.P = _math.sqrt(math.pi * self.R**2 / self.packing_fraction)

        # Areas
        self.area_aperture = math.pi * self.R**2           # cavity opening
        self.area_walls    = 2.0 * math.pi * self.R * self.H  # cylindrical wall
        self.area_base     = math.pi * self.R**2           # flat base

        # dr/dz = 0 (straight cylinder)
        self._dr_dz = 0.0

        # Fix 1 (peer-review): Waveguide modal cutoff for TE11 in a lossy dielectric
        # cylindrical channel.  The PEC approximation λ_c = π·D/j'₁,₁ ≈ 1.706·D
        # assumes infinite conductivity (σ→∞), which is unphysical for alumina
        # (PAA) walls.  The correct cutoff is the complex root of the dielectric
        # boundary characteristic equation solved by solve_te11_mode_complex
        # (waveguide_modes.py, method='perturbation' = first-order lossy correction
        # to the PEC eigenvalue using the actual complex ε_r(ω) of the wall).
        #
        # Reference: Narayanaswamy & Chen, PRB 70, 125101 (2004); Jackson §8.4.
        #
        # The effective cutoff λ_c,eff is the wavelength at which Re(β) → 0.
        # We probe it at a wavelength equal to the PEC estimate (as a seed) and
        # read back the corrected cutoff from the solver.
        _lambda_c_pec = 1.706 * self.diameter_um   # PEC seed (µm)
        try:
            from waveguide_modes import solve_te11_mode_complex
            _modal = solve_te11_mode_complex(
                diameter_um=self.diameter_um,
                wavelength_um=_lambda_c_pec,   # probe near cutoff
                material=self.wall_material,
                method='perturbation',
            )
            # The returned 'cutoff_wavelength_um' is the dielectric-corrected λ_c.
            # For lossy walls the effective cutoff shifts relative to the PEC value.
            _lambda_c_eff = _modal.get('cutoff_wavelength_um', _lambda_c_pec)
            if _lambda_c_eff > 0 and math.isfinite(_lambda_c_eff):
                self.lambda_c_um = float(_lambda_c_eff)
            else:
                self.lambda_c_um = _lambda_c_pec   # fallback
        except Exception:
            # Solver unavailable — fall back to PEC root (backward compatible).
            self.lambda_c_um = _lambda_c_pec
        self.lambda_c = self.lambda_c_um * 1e-6

    def channel_cutoff_wavelength_um(self) -> float:
        """Modal cutoff wavelength (µm). Emission λ > λ_c is evanescent."""
        return self.lambda_c_um

    @property
    def cavity_enhancement(self) -> float:
        """(A_walls + A_base) / A_aperture — should equal alpha_eff / p_esc."""
        return (self.area_walls + self.area_base) / self.area_aperture

    def next_hit(self, pos: np.ndarray, direction: np.ndarray) -> Tuple[float, str, np.ndarray]:
        """Ray intersection: cylindrical cavity interior.

        surface_type in {'aperture', 'wall', 'base'}
        """
        EPS = 1e-13
        t_best = float('inf')
        surf_best = 'none'
        n_x = n_y = n_z = 0.0

        # Aperture at z = H (escapes upward)
        if direction[2] > EPS:
            t_top = (self.H - pos[2]) / direction[2]
            if t_top > EPS and t_top < t_best:
                t_best, surf_best = t_top, 'aperture'
                n_x, n_y, n_z = 0.0, 0.0, 1.0

        # Base at z = 0
        if direction[2] < -EPS:
            t_bot = (0.0 - pos[2]) / direction[2]
            if t_bot > EPS and t_bot < t_best:
                t_best, surf_best = t_bot, 'base'
                n_x, n_y, n_z = 0.0, 0.0, 1.0

        # Cylinder wall: x^2 + y^2 = R^2 (inside looking out)
        dx, dy = direction[0], direction[1]
        px, py = pos[0], pos[1]
        A = dx * dx + dy * dy
        B = 2.0 * (px * dx + py * dy)
        C = px * px + py * py - self.R ** 2
        disc = B * B - 4.0 * A * C

        if A > 1e-30 and disc >= 0.0:
            sq = math.sqrt(disc)
            for sign in (1, -1):   # prefer the smaller positive t
                t_c = (-B + sign * sq) / (2.0 * A)
                if t_c > EPS:
                    z_hit = pos[2] + direction[2] * t_c
                    if 0.0 <= z_hit <= self.H and t_c < t_best:
                        x_hit = px + dx * t_c
                        y_hit = py + dy * t_c
                        t_best, surf_best = t_c, 'wall'
                        # Outward normal for interior hit points TOWARDS cylinder axis (inward)
                        n_x, n_y, n_z = -x_hit / self.R, -y_hit / self.R, 0.0
                    break

        if t_best == float('inf'):
            return float('inf'), 'none', np.zeros(3)
        return t_best, surf_best, np.array([n_x, n_y, n_z])

    def sample_point_on_walls(self) -> Tuple[np.ndarray, np.ndarray]:
        """Uniform random point on cylindrical wall surface."""
        z   = np.random.uniform(0.0, self.H)
        phi = np.random.uniform(0.0, 2 * math.pi)
        x   = self.R * math.cos(phi)
        y   = self.R * math.sin(phi)
        # Normal points inward (toward cavity axis — toward the photon)
        n_out = np.array([-math.cos(phi), -math.sin(phi), 0.0])
        return np.array([x, y, z]), n_out

    def sample_point_on_base(self) -> Tuple[np.ndarray, np.ndarray]:
        """Uniform random point on flat base disk."""
        r   = self.R * math.sqrt(np.random.random())
        phi = np.random.uniform(0.0, 2 * math.pi)
        return (np.array([r * math.cos(phi), r * math.sin(phi), 0.0]),
                np.array([0.0, 0.0, 1.0]))
