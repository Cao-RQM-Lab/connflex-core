# Sidewall capacitance

Where the device layer is thick against its trenches, the field between
facing sidewalls is two-dimensional in the layout plane, and the capacitance
of a layer is that of the plane problem times its thickness. This page solves
that plane problem for a set of conductors by a boundary element method:
charge on straight elements of the outlines, a pair integral in closed form,
the capacitance matrix, its exact first and second derivatives for rigid
motion of the conductors, and the division of a large mesh into windows.
`sidewall_capacitance` runs it on polygons, and the
[capacitance model](capacitance.md) on the electrodes of a structure.

Facing sidewalls of device silicon form capacitors across trenches. Where the
device layer is thick against the width of a trench, the field between
sidewalls does not vary through the thickness, and capacitance is that of a
two-dimensional problem in the layout plane times the thickness. This section
treats the plane problem for a set of conductors: regions bounded by
polygons, conductor $a$ at potential $V_a$, in a medium of permittivity
$\varepsilon = \varepsilon_0 \varepsilon_a$. Every charge and capacitance
here is per unit thickness.

Three methods find the capacitance of the plane problem. The field method
of this section solves the field between the conductors in one dense system
and approaches the exact answer of the plane problem as its mesh is
refined. The same field solved in [windows](#windows) costs time in
proportion to the length of outline instead of its cube. The plates method
of the [next page](methods.md#sidewalls-as-facing-plates) takes every wall that
faces another conductor as a parallel plate with it, costs far less again,
and leaves out all fringe field. Each pair of interacting nets takes the
method that its [feature size](methods.md#method-for-each-pair-of-nets) selects.

Charge sits on the outlines $\Gamma_a$ with density $\sigma$ and produces
the potential

$$
\varphi(\mathbf{r}) = -\frac{1}{2 \pi \varepsilon} \sum_a \oint_{\Gamma_a} \sigma(\mathbf{r}') \ln |\mathbf{r} - \mathbf{r}'|\, ds' + \varphi_0. \tag{61}
$$

In two dimensions the potential of a net charge grows as the logarithm of
distance, so capacitance to a ground at infinity does not exist. The set of
conductors is therefore neutral as a whole, and the constant $\varphi_0$,
the potential far away, is an unknown that floats:

$$
\varphi(\mathbf{r}) = V_a \quad \text{for } \mathbf{r} \text{ on } \Gamma_a, \qquad
\sum_a \oint_{\Gamma_a} \sigma\, ds = 0. \tag{62}
$$

Only differences of potential matter. Equations (61) and (62) have one
solution. Its potential is constant inside every conductor, so $\sigma$ is
the charge density that the field outside leaves on the surface.

### Elements

Outlines are divided into $N$ straight elements. Element $j$ has length
$L_j$ and carries charge $q_j$, spread evenly along it. Requiring the first
condition of (62) on average over each element, which is Galerkin's method
with constant weights, gives

$$
\sum_{j=1}^{N} P_{ij}\, q_j + \varphi_0 = V_{a(i)}, \quad i = 1, \dots, N, \qquad
\sum_{j=1}^{N} q_j = 0, \tag{63}
$$

where $a(i)$ is the conductor of element $i$ and

$$
P_{ij} = -\frac{I_{ij}}{2 \pi \varepsilon\, L_i L_j}, \qquad
I_{ij} = \int_{\text{element } i} \int_{\text{element } j} \ln |\mathbf{r} - \mathbf{r}'|\, ds'\, ds. \tag{64}
$$

$P_{ij}$ is the average potential on element $i$ of a unit charge on element
$j$. It is symmetric, which makes the capacitance matrix symmetric exactly
and not only in the limit of a fine mesh.

### Pair integral

The double integral of (64) has a closed form. Write a point as the complex
number $z = x + \mathrm{i} y$ and element $i$ as $z_i(s) = z_i^0 + \tau_i s$
for $0 \le s \le L_i$, with end points $z_i^0$ and $z_i^1$ and unit tangent
$\tau_i = (z_i^1 - z_i^0) / L_i$. Since $\ln |z| = \operatorname{Re} \ln z$,

$$
I_{ij} = \operatorname{Re} \int_0^{L_i} \int_0^{L_j} \ln \bigl( z_i(s) - z_j(s') \bigr)\, ds'\, ds. \tag{65}
$$

The argument of the logarithm is linear in $s$ and $s'$. With

$$
G(u) = u^2 \left( \tfrac{1}{2} \ln u - \tfrac{3}{4} \right), \qquad G''(u) = \ln u, \tag{66}
$$

$\partial^2 G(z_i(s) - z_j(s')) / \partial s\, \partial s' = -\tau_i \tau_j
\ln(z_i(s) - z_j(s'))$, so the integral is a second difference of $G$ over
the four pairs of end points, $u_{kl} = z_i^k - z_j^l$:

$$
I_{ij} = \operatorname{Re} \left\{ -\frac{1}{\tau_i \tau_j} \bigl[ G(u_{11}) - G(u_{10}) - G(u_{01}) + G(u_{00}) \bigr] \right\}. \tag{67}
$$

The four values $u_{kl}$ are the corners of a parallelogram in the complex
plane, and the logarithm must be one continuous branch over it. Two
different elements that do not cross leave the origin outside the
parallelogram or at one of its corners, so the branch cut can be laid away from it: every $u_{kl}$ is
multiplied by $\bar{u}_c / |u_c|$ before the principal logarithm is taken,
where $u_c$ is the difference of the element centers. The constant this
removes from the logarithm multiplies $L_i L_j$ by a purely imaginary number
and drops out of the real part. Elements that share an end point have
$u_{kl} = 0$ at one corner, where $G(0) = 0$. For an element with itself
(67) reduces to

$$
I_{ii} = L_i^2 \left( \ln L_i - \tfrac{3}{2} \right). \tag{68}
$$

Equation (67) subtracts numbers of size $|u|^2 \ln |u|$ to leave one of size
$L_i L_j \ln |u|$, so it loses accuracy when the elements are short against
their distance. For such pairs the logarithm is expanded about $u_c$
instead. With $\alpha = \tau_i L_i / (2 u_c)$ and
$\beta = \tau_j L_j / (2 u_c)$, odd powers integrate to zero and

$$
I_{ij} = L_i L_j \operatorname{Re} \left[ \ln u_c - \sum_{m=1}^{M} \frac{S_m}{2m} \right], \qquad
S_m = \sum_{k=0}^{m} \binom{2m}{2k} \frac{\alpha^{2k} \beta^{2m-2k}}{(2k+1)(2m-2k+1)}. \tag{69}
$$

The series is used when $L_i + L_j \le 2 \rho\, |u_c|$ with $\rho = 1/8$,
and is cut at $M = 7$. Its remainder is below $\rho^{16} / 16$, about
$2 \times 10^{-16}$ of $L_i L_j$.

### Capacitance matrix

Let $\mathbf{E}$ be the $N \times n$ matrix with $E_{ia} = 1$ when element
$i$ lies on conductor $a$ and zero otherwise, $\mathbf{v}$ the potentials of
the $n$ conductors, and $\mathbf{1}$ the vector of ones. Eliminating
$\varphi_0$ from (63) gives the element charges and the capacitance matrix,

$$
\mathbf{q} = \mathbf{M} \mathbf{E} \mathbf{v}, \qquad
\mathbf{M} = \mathbf{P}^{-1} - \frac{\mathbf{P}^{-1} \mathbf{1} \mathbf{1}^\top \mathbf{P}^{-1}}{\mathbf{1}^\top \mathbf{P}^{-1} \mathbf{1}}, \qquad
\mathbf{c} = \mathbf{E}^\top \mathbf{M} \mathbf{E}, \tag{70}
$$

where $c_{ab}$ is the charge on conductor $a$ per volt on conductor $b$ with
all others at zero. It has three properties by construction:

$$
\mathbf{c}^\top = \mathbf{c}, \qquad \mathbf{c}\, \mathbf{1} = \mathbf{0}, \qquad
\mathbf{v}^\top \mathbf{c}\, \mathbf{v} \ge 0. \tag{71}
$$

The second states that a common potential on all conductors moves no charge,
so $-c_{ab}$ is the capacitance between conductors $a$ and $b$. Stored
energy is $\tfrac{1}{2} \mathbf{v}^\top \mathbf{c}\, \mathbf{v}$. It equals
the largest value of $\mathbf{v}^\top \mathbf{E}^\top \mathbf{q} -
\tfrac{1}{2} \mathbf{q}^\top \mathbf{P} \mathbf{q}$ over neutral charges
$\mathbf{q}$, and charge that is constant on each element is a restricted
choice, so the computed energy never exceeds the exact one for the same
polygons and rises toward it as the mesh is refined. Capacitance of a layer
of thickness $t$ is $t\, \mathbf{c}$.

### Mesh

Charge density is smooth along a straight sidewall, diverges at a convex
corner of a conductor, as $r^{-1/3}$ at a right angle, and vanishes at a
concave one. Elements are therefore at most $h$ long and shorten toward
convex corners. A vertex within $10^{-6} h$ of the next, as rounding leaves
them, is dropped first: an element that short has a row of noise, since its
integrals divide by its length. Let a vertex turn the outline by the angle $\theta$ toward
the conductor, so that $\theta > 0$ at a convex corner and $\theta = \pi/2$
at a right angle. The length asked for at the vertex is

$$
h_c = h\, g^{-w}, \qquad w = \min \bigl( \max( 2 \theta / \pi,\, 0 ),\, 1 \bigr), \tag{72}
$$

which is $h / g$ at a right angle or sharper, and $h$ at a straight or
concave vertex, so the many vertices of a drawn arc cause no refinement.
Along an edge of length $\ell$ between vertices that ask for $h_L$ and
$h_R$, the length asked for at distance $x$ from the first grows with
distance from either end,

$$
s(x) = \min ( h,\; h_L + x,\; h_R + \ell - x ). \tag{73}
$$

The edge gets $n$ elements with nodes $x_k$ at equal steps of the element
count that $s$ implies,

$$
\nu(x) = \int_0^x \frac{dx'}{s(x')}, \qquad n = \lceil \nu(\ell) \rceil, \qquad
\nu(x_k) = \frac{k}{n}\, \nu(\ell), \quad k = 0, \dots, n. \tag{74}
$$

On a rising ramp $\nu = \ln(1 + x / h_L)$, so lengths there grow by a factor
of at most $e$ from one element to the next, and no element is longer than
$h$. An edge shorter than the lengths asked for at its ends is one element.

### Rigid motion

Conductors are grouped into rigid bodies. Body $B$ travels by $(a_x, a_y)$
and turns by $a_\psi$ about a point $\mathbf{c}_B$, written $c_B$ as a
complex number, which to first order moves a point $z$ of its outlines by

$$
\delta z = a_x + \mathrm{i}\, a_y + \mathrm{i}\, a_\psi\, (z - c_B). \tag{75}
$$

Elements keep their lengths, so $\mathbf{E}$ and $\mathbf{1}$ in (70) do not
change and the capacitance matrix depends on the motion through
$\mathbf{P}$ alone. $\mathbf{M}$ is the leading block of the inverse of the
bordered matrix of (63), whose border is constant, so for any parameter $p$
of the motion

$$
\frac{\partial \mathbf{c}}{\partial p} = -\mathbf{X}^\top \frac{\partial \mathbf{P}}{\partial p}\, \mathbf{X}, \qquad
\mathbf{X} = \mathbf{M} \mathbf{E}, \tag{76}
$$

where column $a$ of $\mathbf{X}$ holds the element charges for one volt on
conductor $a$. No further system is solved. The pair integral (65) depends
only on how two elements lie against each other, so its derivative is zero
when both belong to one body. For element $i$ on the moving body and
element $j$ on another, differentiating (65) under the integral with (75)
brings in

$$
T_{ij} = \int_0^{L_i} \int_0^{L_j} \frac{ds'\, ds}{z_i(s) - z_j(s')}, \qquad
W_{ij} = \int_0^{L_i} \int_0^{L_j} \frac{z_i(s) - c_B}{z_i(s) - z_j(s')}\, ds'\, ds, \tag{77}
$$

and gives

$$
\frac{\partial I_{ij}}{\partial a_x} = \operatorname{Re} T_{ij}, \qquad
\frac{\partial I_{ij}}{\partial a_y} = -\operatorname{Im} T_{ij}, \qquad
\frac{\partial I_{ij}}{\partial a_\psi} = -\operatorname{Im} W_{ij}. \tag{78}
$$

Both integrals are second differences like (67). With
$H(u) = G'(u) = u\, (\ln u - 1)$, the same branch of the logarithm, signs
$\epsilon_{kl} = +1$ for $k = l$ and $-1$ otherwise, and $I^\mathrm{c}_{ij}$
the complex number in braces in (67),

$$
T_{ij} = -\frac{1}{\tau_i \tau_j} \sum_{k, l} \epsilon_{kl}\, H(u_{kl}), \qquad
W_{ij} = -\frac{1}{\tau_i \tau_j} \sum_{k, l} \epsilon_{kl}\, H(u_{kl})\, (z_i^k - c_B) - I^\mathrm{c}_{ij}. \tag{79}
$$

For separated pairs the expansion that gave (69) gives, with $z_i^c$ the
center of element $i$,

$$
T_{ij} = \frac{L_i L_j}{u_c} \left[ 1 + \sum_{m=1}^{M} S_m \right], \qquad
W_{ij} = (z_i^c - c_B)\, T_{ij} - L_i L_j \sum_{m=1}^{M} R_m, \qquad
R_m = \sum_{k=1}^{m} \binom{2m-1}{2k-1} \frac{\alpha^{2k} \beta^{2m-2k}}{(2k+1)(2m-2k+1)}, \tag{80}
$$

in which the first term of $W_{ij}$ is the pull on the element acting at
its center and the second its turning about that center. Since
$I_{ij} = I_{ji}$, a pair whose second element moves is covered by swapping
the two. With $s_{ia} = X_{ia} / L_i$ the charge density on element $i$ per
volt on conductor $a$, (64) and (76) combine to

$$
\frac{\partial c_{ab}}{\partial p_B} = \frac{1}{2 \pi \varepsilon} \sum_{i \in B} \sum_{j \notin B}
\left( s_{ia}\, s_{jb} + s_{ib}\, s_{ja} \right) \frac{\partial I_{ij}}{\partial p_B} \tag{81}
$$

for each of the three parameters $p_B$ of body $B$.

Equation (81) is the derivative of the computed matrix (70) itself, not a
separate approximation of the derivative of the exact one, so force and
energy belong to one function of position: with potentials $\mathbf{v}$ the
force on body $B$ is $\tfrac{1}{2} \mathbf{v}^\top (\partial \mathbf{c} /
\partial p_B)\, \mathbf{v}$ per unit thickness, and it equals the change of
$\tfrac{1}{2} \mathbf{v}^\top \mathbf{c}\, \mathbf{v}$ when the bodies and
their elements are moved, to rounding. It is the pull of the charge on other
bodies on the charge of body $B$. The field of other bodies is smooth along
the outlines of $B$, so the squares of charge density that diverge at
corners, which a surface stress integral would need, never appear. Moving
every body together changes nothing:

$$
\sum_B \frac{\partial \mathbf{c}}{\partial a_x^B} = \sum_B \frac{\partial \mathbf{c}}{\partial a_y^B} = \mathbf{0}, \qquad
\sum_B \frac{\partial \mathbf{c}}{\partial a_\psi^B} = \mathbf{0} \quad \text{for a common } \mathbf{c}_B, \tag{82}
$$

which holds pair by pair, because $T_{ij} + T_{ji} = 0$ and
$W_{ij} + W_{ji} = L_i L_j$ is real: forces and torques on the bodies
balance.

### Second derivatives

Electrostatic stiffness needs the second derivative of capacitance, for
which (75) must be the rigid motion it approximates. Body $B$ travels by
$a = a_x + \mathrm{i}\, a_y$ and turns by $a_\psi$ about a point that travels
with it:

$$
z \mapsto c_B + a + e^{\mathrm{i} a_\psi} (z - c_B). \tag{83}
$$

Differentiating (76) once more, with $\partial \mathbf{M} / \partial q =
-\mathbf{M}\, (\partial \mathbf{P} / \partial q)\, \mathbf{M}$, gives

$$
\frac{\partial^2 \mathbf{c}}{\partial p\, \partial q}
= -\mathbf{X}^\top \frac{\partial^2 \mathbf{P}}{\partial p\, \partial q} \mathbf{X}
+ \mathbf{Y}_p^\top \mathbf{M} \mathbf{Y}_q + \mathbf{Y}_q^\top \mathbf{M} \mathbf{Y}_p, \qquad
\mathbf{Y}_p = \frac{\partial \mathbf{P}}{\partial p} \mathbf{X}. \tag{84}
$$

The first term is the change of the pull between fixed charges. The other
two are the pull that arises because the charges themselves move along the
outlines when a body moves; they need $\mathbf{M} \mathbf{Y}_q$, one more
pair of triangular solves per parameter with the factor already at hand.

Take a pair with element $i$ on body $B$ and element $j$ on body $B'$, and
the six parameters $(a_x, a_y, a_\psi)$ of $B$ and $(a_x', a_y', a_\psi')$
of $B'$. With $r_k = z_i^k - c_B$ and $r_l' = z_j^l - c_{B'}$, the end-point
differences $u_{kl}$ of (67) move under (83) with

$$
\frac{\partial u_{kl}}{\partial p} = v_p = \bigl( 1,\; \mathrm{i},\; \mathrm{i}\, r_k,\; -1,\; -\mathrm{i},\; -\mathrm{i}\, r_l' \bigr), \qquad
\frac{\partial^2 u_{kl}}{\partial a_\psi^2} = -r_k, \qquad
\frac{\partial^2 u_{kl}}{\partial a_\psi'^2} = r_l', \tag{85}
$$

all other second derivatives being zero, and the product $\tau_i \tau_j$
turns by $e^{\mathrm{i} (a_\psi + a_\psi')}$. Let $\varrho_p$ be one for the
two rotations and zero for the four translations. Differentiating the
complex integral $I^\mathrm{c}_{ij}$ of (67) gives

$$
\begin{aligned}
\frac{\partial I^\mathrm{c}_{ij}}{\partial p} &= \Phi_p - \mathrm{i}\, \varrho_p\, I^\mathrm{c}_{ij}, \qquad
\Phi_p = -\frac{1}{\tau_i \tau_j} \sum_{k,l} \epsilon_{kl}\, H(u_{kl})\, v_p, \\
\frac{\partial^2 I^\mathrm{c}_{ij}}{\partial p\, \partial q} &= -\frac{1}{\tau_i \tau_j} \sum_{k,l} \epsilon_{kl}
\left[ \ln (u_{kl})\, v_p v_q + H(u_{kl}) \frac{\partial^2 u_{kl}}{\partial p\, \partial q} \right]
- \mathrm{i}\, \varrho_p \Phi_q - \mathrm{i}\, \varrho_q \Phi_p - \varrho_p \varrho_q\, I^\mathrm{c}_{ij},
\end{aligned} \tag{86}
$$

whose real parts are the derivatives of $I_{ij}$. The first line contains
(78) and (79): $\Phi_{a_x} = T_{ij}$ and
$\Phi_{a_\psi} - \mathrm{i} I^\mathrm{c}_{ij} = \mathrm{i}\, W_{ij}$.

For separated pairs write (69) as $I^\mathrm{c}_{ij} = L_i L_j\, F$ with
$F = \ln u_c - \sum_{m,k} e_{mk}\, t_{mk} / (2m)$, where
$t_{mk} = \alpha^{2k} \beta^{2m-2k}$ and $e_{mk}$ is its coefficient in
$S_m$. $F$ depends on the motion through $u_c$ and through the two angles,
since $\alpha^2$ turns by $e^{2 \mathrm{i} a_\psi}$ and $\beta^2$ by
$e^{2 \mathrm{i} a_\psi'}$. Its partial derivatives are sums over the same
terms,

$$
\begin{aligned}
F_u &= \frac{1}{u_c} \Bigl[ 1 + \sum e_{mk} t_{mk} \Bigr], &
F_{uu} &= -\frac{1}{u_c^2} \Bigl[ 1 + \sum (2m+1)\, e_{mk} t_{mk} \Bigr], \\
F_\psi &= -\mathrm{i} \sum \frac{k}{m}\, e_{mk} t_{mk}, &
F_{\psi'} &= -\mathrm{i} \sum \frac{m-k}{m}\, e_{mk} t_{mk}, \\
F_{u\psi} &= \frac{2 \mathrm{i}}{u_c} \sum k\, e_{mk} t_{mk}, &
F_{u\psi'} &= \frac{2 \mathrm{i}}{u_c} \sum (m-k)\, e_{mk} t_{mk}, \\
F_{\psi\psi} &= \sum \frac{2 k^2}{m}\, e_{mk} t_{mk}, &
F_{\psi'\psi'} &= \sum \frac{2 (m-k)^2}{m}\, e_{mk} t_{mk}, \qquad
F_{\psi\psi'} = \sum \frac{2 k (m-k)}{m}\, e_{mk} t_{mk}.
\end{aligned} \tag{87}
$$

The difference of centers $u_c$ moves as in (85) with the element centers in
place of the end points; call its first derivatives $w_p$. Writing
$\omega_p$ for the angle that parameter $p$ is, if it is one,

$$
\frac{\partial F}{\partial p} = F_u\, w_p + F_{\omega_p}, \qquad
\frac{\partial^2 F}{\partial p\, \partial q} = F_{uu}\, w_p w_q + F_u \frac{\partial^2 u_c}{\partial p\, \partial q}
+ F_{u \omega_p} w_q + F_{u \omega_q} w_p + F_{\omega_p \omega_q}, \tag{88}
$$

where a term with $\omega_p$ is absent when $p$ is a translation. The first
of these contains (80).

With these, the first term of (84) is a sum over pairs of elements on
different bodies, as (81) is,

$$
-\Bigl[ \mathbf{X}^\top \frac{\partial^2 \mathbf{P}}{\partial p\, \partial q} \mathbf{X} \Bigr]_{ab}
= \frac{1}{2 \pi \varepsilon} \sum_{B < B'} \sum_{i \in B} \sum_{j \in B'}
\left( s_{ia}\, s_{jb} + s_{ib}\, s_{ja} \right) \frac{\partial^2 I_{ij}}{\partial p\, \partial q}, \tag{89}
$$

in which a pair contributes when $p$ and $q$ are each a parameter of one of
its two bodies, and

$$
\bigl[ \mathbf{Y}_p \bigr]_{ia} = -\frac{1}{2 \pi \varepsilon\, L_i} \sum_j \frac{\partial I_{ij}}{\partial p}\, s_{ja}, \tag{90}
$$

summed over the elements $j$ that do not share a body with element $i$.
Translating every body together leaves the gradient unchanged, so the
second derivatives for travel along one axis sum to zero over the bodies of
either parameter.

### Windows

One dense solution of (70) is out of reach for a comb of thousands of
fingers. Charge on a sidewall, though, is set by its surroundings: across a
narrow gap by the wall it faces, and hardly at all by conductors far off. A
mesh with more elements than `max_elements` is therefore solved piece by
piece.

The elements are divided into cores, every element in exactly one, by
halving the set at the median of the element centers along the longer side
of their bounding box until no core holds more than $n_w$ elements. Let
$K_k$ be the elements of core $k$ and $B_k$ the bounding box of their
centers $\mathbf{r}_i$. Window $k$ holds the core and every element within
the overlap $\delta_k$ of that box,

$$
W_k = \bigl\{ i : \mathbf{r}_i \in B_k \oplus [-\delta_k, \delta_k]^2 \bigr\} \supseteq K_k, \tag{91}
$$

and is solved as a problem (63) of its own: its elements alone, each at the
potential of its conductor, neutral together, with a floating potential of
its own. Let $\mathbf{E}_k$ and $\mathbf{M}_k$ be the matrices of (70) for
the elements of $W_k$, and $\mathbf{E}_k^c$ equal to $\mathbf{E}_k$ with the
rows of elements outside the core set to zero. The window supplies the
charge on its core alone,

$$
\mathbf{c}^{(k)} = \mathbf{E}_k^{c\,\top} \mathbf{M}_k \mathbf{E}_k, \tag{92}
$$

whose entry $(a, b)$ is the charge on the core elements of conductor $a$ per
volt on conductor $b$. Every element is in one core, so the sum over windows
counts the charge on every element once:

$$
\hat{\mathbf{c}} = \sum_k \mathbf{c}^{(k)}, \qquad \hat{\mathbf{c}}\, \mathbf{1} = \mathbf{0}. \tag{93}
$$

Rows sum to zero because $\mathbf{M}_k \mathbf{1} = \mathbf{0}$ in every
window. $\hat{\mathbf{c}}$ is not symmetric, since each window sees its own
part of the field, and the capacitance matrix is

$$
c_{ab} = \frac{\hat{c}_{ab} + \hat{c}_{ba}}{2} \quad (a \ne b), \qquad
c_{aa} = -\sum_{b \ne a} c_{ab}, \tag{94}
$$

which has the first two properties of (71). The asymmetry
$\max |\hat{c}_{ab} - \hat{c}_{ba}| / \max |\hat{c}_{ab}|$ is zero when one
window holds the whole mesh. Otherwise it measures what the windows leave
out, and it is reported with the result.

Elements keep their cores as bodies move, and a window its elements while
(164) holds, so (92) depends on the motion through $\mathbf{P}_k$ alone, as
(70) does. With
$\mathbf{X}_k = \mathbf{M}_k \mathbf{E}_k$ and
$\mathbf{Y}_k = \mathbf{M}_k \mathbf{E}_k^c$,

$$
\begin{aligned}
\frac{\partial \mathbf{c}^{(k)}}{\partial p} &= -\mathbf{Y}_k^\top \frac{\partial \mathbf{P}_k}{\partial p}\, \mathbf{X}_k, \\
\frac{\partial^2 \mathbf{c}^{(k)}}{\partial p\, \partial q} &= -\mathbf{Y}_k^\top \frac{\partial^2 \mathbf{P}_k}{\partial p\, \partial q}\, \mathbf{X}_k
+ \Bigl( \frac{\partial \mathbf{P}_k}{\partial p} \mathbf{Y}_k \Bigr)^{\!\top} \mathbf{M}_k\, \frac{\partial \mathbf{P}_k}{\partial q} \mathbf{X}_k
+ \Bigl( \frac{\partial \mathbf{P}_k}{\partial q} \mathbf{Y}_k \Bigr)^{\!\top} \mathbf{M}_k\, \frac{\partial \mathbf{P}_k}{\partial p} \mathbf{X}_k,
\end{aligned} \tag{95}
$$

which are (76) and (84) with $\mathbf{Y}_k$ in place of $\mathbf{X}_k$ on
the left, evaluated with the same pair derivatives (85) to (88). They pass
through (93) and (94) as $\mathbf{c}^{(k)}$ does and are the exact
derivatives of the matrix (94).

The overlap is found for each window. It starts at $\delta_k = \min(4 h, d_i)$ and
doubles until the charge on the core has settled,

$$
\max_{a,b} \bigl| c^{(k)}_{ab}(\delta_k) - c^{(k)}_{ab}(\delta_k / 2) \bigr| < \tau\, \max_{a,b} \bigl| c^{(k)}_{ab}(\delta_k) \bigr|, \tag{96}
$$

with $\tau$ the window tolerance, or until $\delta_k$ reaches the interaction
distance $d_i$ when one is given, or until the next window would hold more
than `max_elements`. A window whose elements all lie on one conductor
supplies nothing and is not solved; its overlap still doubles, since a
second conductor may come into view.

Measured against one dense solution of the same mesh, a comb of 30 fingers
with gaps of 3 um, 4825 elements in 32 windows:

| $d_i$ | Capacitance | Force along fingers | Stiffness across fingers | Asymmetry |
| --- | --- | --- | --- | --- |
| 30 um | $-2.6 \times 10^{-4}$ | $-6.6 \times 10^{-3}$ | $-4.1 \times 10^{-5}$ | $4.1 \times 10^{-4}$ |
| 60 um | $+3.2 \times 10^{-4}$ | $-7.4 \times 10^{-4}$ | $+2.3 \times 10^{-6}$ | $5.2 \times 10^{-5}$ |

Charge in the gaps settles within an overlap of two gaps. Force along the
fingers converges slowest: it comes from the charge at finger tips, which
face open space, where the plane field reaches far. $d_i$ belongs at a
few times the widest trench that should count: a window whose overlap
stops short of two trench widths leaves out part of the field across it.

!!! note "Scaling and factorization"
    Lengths are divided by twice the diagonal of the bounding box of all
    outlines before (64) is formed. Every distance is then below one, which
    makes $\mathbf{P}$ positive definite, and (70) is solved with one
    Cholesky factorization and $n + 1$ right-hand sides. The scaling adds
    the same constant to every $P_{ij}$, which cancels through the second
    condition of (63), and $\varepsilon$ is applied to the result.

!!! note "Holes that hold nothing"
    A hole in a conductor is meshed only when another polygon lies inside
    it. A cavity that holds no conductor holds no field and its wall carries
    no charge.

!!! note "Conductors must not touch"
    Conductors must not overlap or touch: two regions in contact at
    different potentials have no finite capacitance. Touching regions at
    one potential are passed as one conductor.

!!! note "Accuracy of the pair integral"
    In (67) the remaining loss of accuracy is set by how short an element is
    against a neighbor, which (72) bounds by $g$. Measured against 40-digit
    arithmetic over random pairs with length ratios up to 64, both forms
    agree with the exact integral to $10^{-12}\, L_i L_j$.

!!! note "Cost of one dense solution and of windows"
    One solution is dense: memory grows as $N^2$ and time as $N^3$, and
    gradient and second derivatives add one pass over the pairs of elements
    that lie on different bodies, about as long as forming $\mathbf{P}$. In
    windows the cost grows with $N$ alone and memory with the size of one
    window, so a cluster far too large to hold as one system is solved in
    the time of its windows. Every window is scaled by its own bounding box.

!!! note "Why a window floats"
    A window is neutral and floats, as the whole problem does, because no
    other condition at its cut edge is right in the plane: charge in a gap
    depends only on the difference of the potentials across it, which a
    floating window keeps, while charge on a wall that faces open space
    depends on the potential far away, which a window can only guess. That
    is the part $d_i$ and $\tau$ trade against time. When $d_i$ is given the
    overlap stops there, the distance beyond which the plane field is not
    trusted anyway.

!!! note "A first window too large"
    A first window that already holds more than `max_elements` is an error.

!!! note "Cores beside rows"
    Where rows of identical cells, below, have taken elements out of the
    cores, the elements left are scattered, and a core is also halved while
    its first window would hold more than twice $n_w$, so that no window
    grows to span the rows between its pieces.

!!! note "Windows that are translates"
    A window that is a translate of one solved before takes its solution:
    element for element, with the same elements in its core and its
    conductors and bodies in the same places, it holds the same field.
    Every window finds its derivatives about the lowest corner of the box
    $B_k$ of its core and moves them to the pivots of its bodies with the
    lever arms of (155), so that the copies differ only by those arms.
    Translates are found by the coordinates of the elements relative to
    that corner, rounded to $10^{-11}$ m, and cores are cut on centers
    rounded the same way, so that groups of elements that are translates
    are cut into cores that are. Windows cut by medians over a whole mesh
    are rarely translates; those of the groups of (156) often are.

!!! note "Processes, threads, and compiled kernels"
    Windows and the cross-sections of the fringe are independent pieces of
    work, and worker processes solve their fields: the numpy work around
    the kernels holds the interpreter's lock, which threads cannot share.
    Threads search the windows, cut the fringe lines across the outlines,
    and assemble the blocks of rows of one dense system, in shapely and in
    compiled code that release the lock. A translate met by two windows at
    once is solved by the first and taken by the second, and the two differ
    by rounding of their coordinates. Sums over windows are taken in the
    order of their cores, so the number of workers changes results only
    within that rounding. The pair integrals (67) and (69) and their
    derivatives (85) to (88) are evaluated one pair at a time in code
    compiled by numba when it is installed; otherwise in numpy, in passes
    over blocks of pairs.

!!! note "The pivot travels with the body"
    The point a body turns about travels with the body, as a center of mass
    does. A second derivative for travel and rotation taken about a point
    fixed in space differs from this one by the gradient.

!!! note "Which conductors enter"
    A conductor couples to every other conductor in the plane, however far,
    because the plane field decays only as the logarithm of distance. Which
    conductors enter one solution is the caller's choice.

## Symbols on this page

| Symbol | Meaning | In the code |
| --- | --- | --- |
| $\varepsilon$ | permittivity of the medium in the trenches | `ElectrostaticParams.ambient_permittivity` times `electrostatics.VACUUM_PERMITTIVITY` |
| $h$ | longest element | `ElectrostaticParams.element_size` |
| $g$ | refinement of elements at convex corners | `ElectrostaticParams.corner_grading` |
| $L_j$, $q_j$ | length and charge of element $j$ | `SidewallCapacitance.mesh`, `SidewallCapacitance.charge` |
| $\mathbf{P}$, $I_{ij}$ | potential matrix and pair integral | `electrostatics.dense.solve_charges`, `electrostatics.kernel.pair_integrals` |
| $\rho$, $M$ | reach and order of the far-field series of the pair integral | `electrostatics.kernel.SERIES_REACH`, `SERIES_ORDER` |
| $\mathbf{c}$ | capacitance matrix per unit thickness | `SidewallCapacitance.capacitance` |
| $a_x$, $a_y$, $a_\psi$ | travel and turn of a body in the plane | last axes of `SidewallCapacitance.gradient` and `hessian` |
| $c_B$ | point a body turns about | `SidewallCapacitance.centers` |
| $T_{ij}$, $W_{ij}$ | derivatives of the pair integral for travel and turn | `electrostatics.kernel.pair_derivatives` |
| $n_w$ | most elements in the core of a window | `ElectrostaticParams.window_elements` |
| $\delta_k$, $\tau$ | overlap of a window and the tolerance that stops its growth | `electrostatics.windows.FIRST_OVERLAP`, `ElectrostaticParams.window_tolerance` |
| $\hat{\mathbf{c}}$ and its asymmetry | sum over windows before symmetrization, and the mismatch it leaves | `SidewallCapacitance.asymmetry`, `SidewallCapacitance.windows` |
| $d_i$ | interaction distance, the widest overlap of a window | `ElectrostaticParams.interaction_distance` |
