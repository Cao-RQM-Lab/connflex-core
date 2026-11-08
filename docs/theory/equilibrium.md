# Static equilibrium

With potentials on its nets, a structure comes to rest where its beams balance
the electrostatic forces. This page poses that balance, solves it by Newton
iteration on the capacitance model evaluated at every iterate, and gives the
measures that go with a solution: whether it is stable, how far it is from
pull-in, how close its gaps are to touching, and how a net that holds a fixed
charge enters. `Structure.equilibrium` solves one operating point and returns
an `Equilibrium`.

With potentials on its nets a structure comes to rest where the pull of its
beams balances the electrostatic forces. Nothing flows at rest, so no link
carries a voltage and the electrodes of a net share its potential:
$\mathbf{v} = \mathbf{S}\, \mathbf{v}^\text{net}$ with $\mathbf{S}$ of (113).
Beams are linear, with the stiffness $\mathbf{K}^\text{full}$ of (35), and
the electrostatic force is that of (119) at the displaced bodies. With
mechanical loads $\mathbf{f}_0$ on the bodies, the displacement $\mathbf{q}$
solves

$$
\mathbf{K}^\text{full} \mathbf{q} = \mathbf{f}(\mathbf{q}) + \mathbf{f}_0, \qquad
f_k(\mathbf{q}) = \frac{1}{2}\, \mathbf{v}^\top \frac{\partial \mathbf{C}}{\partial q_k}(\mathbf{q})\, \mathbf{v}. \tag{120}
$$

The force grows as gaps close, so (120) is not linear and is solved by
Newton iteration from a start $\mathbf{q}_0$, rest unless another is given:

$$
\bigl( \mathbf{K}^\text{full} + \mathbf{K}_e(\mathbf{q}_i) \bigr)\, \Delta \mathbf{q}
= \mathbf{f}(\mathbf{q}_i) + \mathbf{f}_0 - \mathbf{K}^\text{full} \mathbf{q}_i, \qquad
\mathbf{q}_{i+1} = \mathbf{q}_i + \Delta \mathbf{q}, \tag{121}
$$

with $\mathbf{K}_e$ of (119). Every iterate evaluates the capacitance model
at $\mathbf{q}_i$ for its gradient and second derivative, which are exact
there, so the iteration converges quadratically. Displacements mix meters
and radians, and steps are therefore measured by the strain energy they
would store in the beams. Iteration stops when

$$
\Delta \mathbf{q}^\top \mathbf{K}^\text{full} \Delta \mathbf{q}
\le \tau_s^2\; \mathbf{q}_{i+1}^\top \mathbf{K}^\text{full} \mathbf{q}_{i+1}, \tag{122}
$$

with $\tau_s$ the tolerance of the solve, or one step after the left side
has fallen below $\tau_s$ times the right side without the square. A step
that small is followed by one of about its square, which meets (122), unless
rounding in the forces keeps it larger, and then further steps would not
shrink it either.

Equation (120) has more than one solution in general. Below pull-in a plate
on a spring has two, one that holds and one that tips over; a toothed rotor
has a detent on either side of a balance point. Newton iteration converges
to the solution near its start, so the start chooses among them. A device
reaches its state by its potentials changing from earlier ones, and the
equilibrium it sits in is the one that grows continuously out of the
equilibrium it had. Each solve finds one equilibrium at one set of
potentials; a sequence, such as potentials raised from zero or phases of a
motor switched in turn, is a series of solves, each started from the
equilibrium of the one before at nearby potentials.

An equilibrium holds only if a small push meets a restoring force in every
direction. The tangent stiffness $\mathbf{K}^\text{full} +
\mathbf{K}_e(\mathbf{q})$ answers that, and against the mass matrix of (34)
its smallest eigenvalue,

$$
\bigl( \mathbf{K}^\text{full} + \mathbf{K}_e(\mathbf{q}) \bigr)\, \boldsymbol{\phi}
= \lambda\, \mathbf{M}^\text{full} \boldsymbol{\phi}, \qquad
\lambda_\text{min} > 0, \tag{123}
$$

is the squared angular frequency of the softest motion under these
potentials, as in (38). The equilibrium is stable while it is positive.
Against the beam stiffness instead,

$$
\bigl( \mathbf{K}^\text{full} + \mathbf{K}_e(\mathbf{q}) \bigr)\, \boldsymbol{\psi}
= \mu\, \mathbf{K}^\text{full} \boldsymbol{\psi}, \qquad
\mu_\text{min} = \min_{\boldsymbol{\psi}} \frac{\boldsymbol{\psi}^\top \bigl( \mathbf{K}^\text{full} + \mathbf{K}_e \bigr) \boldsymbol{\psi}}{\boldsymbol{\psi}^\top \mathbf{K}^\text{full} \boldsymbol{\psi}}, \tag{124}
$$

the smallest eigenvalue is the share of beam stiffness left in the
direction the field softens most: 1 without potentials, 0 at pull-in, and
of the same sign as $\lambda_\text{min}$. The two answer different
questions. The motion that the field softens most can be stiff to begin
with, such as a plate moving toward the handle, and its frequency stays
far above that of the softest motion until close to pull-in, while
$\mu_\text{min}$ falls steadily. A square plate on four beams over the
handle at 95 percent of its pull-in voltage has $\mu_\text{min} = 0.44$,
while its softest motion, a travel in the plane that the field does not
touch, keeps the frequency it has without voltage.
For a plate of area $A$ on a spring $k$ over a gap $g_0$, (120) reads
$k\, z = -\varepsilon A V^2 / \bigl( 2 (g_0 + z)^2 \bigr)$ and (123) is
$k > \varepsilon A V^2 / (g_0 + z)^3$.

### Pull-in

Raised from zero, the potentials carry the equilibrium that grows out of
rest until it ends, in one of two ways. At a fold $\lambda_\text{min}$
falls to zero and no equilibrium exists nearby beyond it: a plate on a
spring meets its fold at one third of its gap. At a loss of symmetry the
equilibrium goes on but no longer holds against a push: a comb between two
stators at one potential stays centered while the spring across its fingers
is used up by (117). The two closed forms are

$$
V^* = \sqrt{\frac{8\, k\, g_0^3}{27\, \varepsilon A}} \ \text{ at } z = -\frac{g_0}{3}, \qquad
V^*_\perp = \sqrt{\frac{k_\perp\, g^3}{N \varepsilon\, t_\text{device}\, \ell}}. \tag{125}
$$

Near a fold $\mu_\text{min}$ falls toward zero, and beyond it Newton
iteration finds no solution near its start. At a loss of symmetry the
solution it finds has $\mu_\text{min} < 0$ and is reported as not stable.
Locating pull-in is therefore a search over potentials by repeated solves,
which a caller writes; the usage page shows one.

### Contact

Nothing in the model describes two bodies that touch, so an equilibrium is
meaningful only while every gap between two bodies stays open. The clearance
of a displacement is the smallest ratio of such a gap to its width at rest,

$$
\chi(\mathbf{q}) = \min \Bigl\{
\min_{\text{plates}}\, \min_{\mathbf{r} \in V_{ab}} \frac{g_0 + h_0 + h_x\, r_x + h_y\, r_y}{g_0},\;
\min_{(a, b)} \frac{\operatorname{dist}\bigl( E_a(\mathbf{q}), E_b(\mathbf{q}) \bigr)}{\operatorname{dist}(E_a, E_b)}
\Bigr\}. \tag{126}
$$

The first runs over plates between two bodies, with the rise $h_0$ and
slopes $h_x$, $h_y$ of (45) and $V_{ab}$ the vertices of the convex hull of
the overlap $G_{ab}$, measured from the reference point: the gap (45) is a
plane, so over the overlap it is smallest at one of them, and its numerator
is the gap at that corner. The second runs over pairs of device electrodes
that move with different bodies and lie within $d_i$ of each other at rest
or at $\mathbf{q}$, with the outlines carried by (118), and is the exact
distance between their silicon; a pair that was not within $d_i$ at rest is
in no cluster, so its capacitance is not in the model, but its gap is.
Electrodes that touch are an error. The model at every
equilibrium carries $\chi$ and the two electrodes across that gap; how
small a clearance to accept is the caller's judgment.

### Floating nets

A net that nothing drives keeps the charge it has. Split the nets into
driven ones, $d$, whose potentials are given, and floating ones, $f$, whose
charges $\mathbf{Q}_f$ are given. With the matrix between nets
$\mathbf{C}^\text{net}(\mathbf{q})$ of (113) in blocks, the charge on the
floating nets fixes their potentials at every displacement:

$$
\mathbf{C}^\text{net}_{ff}\, \mathbf{v}_f + \mathbf{C}^\text{net}_{fd}\, \mathbf{v}_d = \mathbf{Q}_f
\quad \Longrightarrow \quad
\mathbf{v}_f(\mathbf{q}) = \mathbf{C}_{ff}^{\text{net}\,-1} \bigl( \mathbf{Q}_f - \mathbf{C}^\text{net}_{fd}\, \mathbf{v}_d \bigr). \tag{127}
$$

The force keeps the form (120), now with these potentials. Its change with
displacement gains a part, because moving a body changes
$\mathbf{v}_f$: with $\boldsymbol{\Gamma}_f$ the rows for floating nets of
$\Gamma_{\nu k} = \sum_\mu \bigl( \partial C^\text{net}_{\nu\mu} / \partial q_k \bigr) v_\mu$,

$$
\mathbf{K}_t = \mathbf{K}^\text{full} + \mathbf{K}_e
+ \boldsymbol{\Gamma}_f^\top\, \mathbf{C}_{ff}^{\text{net}\,-1}\, \boldsymbol{\Gamma}_f \tag{128}
$$

replaces the tangent stiffness in (121), (123), and (124). The added part is never
negative: a body that moves toward a floating electrode lowers its
potential, which weakens the pull. A plate that holds a fixed charge and
has no other capacitance loses no stiffness at all, and never pulls in.
$\mathbf{C}^\text{net}_{ff}$ must be positive definite, which it is when
every floating net has capacitance to a driven one.

!!! note "Halving a step into a gap"
    A Newton step that would close a gap between a plate and the handle,
    or bring electrodes of two bodies into contact, is halved until it does
    not, up to twelve times. A halved step never ends the iteration with an
    equilibrium.

!!! note "When the iteration gives up"
    An iteration that has not ended after `max_iterations`, whose steps
    were halved or did not shrink `patience` times in a row, eight by
    default, whose step closed a
    gap through all twelve halvings, whose tangent stiffness is singular, or
    that meets a force or stiffness that is not finite, as a model evaluated
    far beyond its range can return, finds no
    equilibrium. Beyond a fold that is the truth, and giving up early saves
    most of the cost of the attempt. It also happens with the plates method
    on teeth that stand tip to tip: the force jumps at the cusp of (97) and
    the iteration steps back and forth across it at any potential. Such a
    device needs the field method.

!!! note "A plate that closes at a corner"
    A plate whose gap (45) closes at a corner of its overlap is an error
    when evaluated, as electrodes of one cluster that touch are, so a Newton
    step that would close it is halved like one into contact. A tilted plate
    is also where the capacitance of (55), of second order in tilt, is least
    accurate: at a clearance of one half its integrand is 12 percent low at
    the closest corner.

!!! note "Bounds spare most distances"
    A body moved by $(a_x, a_y, a_\psi)$ of (118) carries no point of an
    electrode farther than $\sqrt{a_x^2 + a_y^2} + |a_\psi|\, R$, with $R$
    the distance from the body's center to the farthest corner of the
    electrode's bounding box, so a gap closes by at most the sum of the two
    electrodes' bounds. A pair whose gap cannot have fallen below the
    smallest ratio found so far is not measured, a cluster whose pairs
    cannot have closed is not checked for electrodes that touch, and of the
    pairs that no cluster measures only those that may have come within the
    interaction distance are carried and measured. Those pairs are found
    once at rest, within the interaction distance plus twice the largest
    bound, and found again at least twice as far out when a larger motion
    could bring farther pairs within reach. The values are those of the
    full check.

!!! note "A failure from far away"
    An iteration can fail from a distant start where it succeeds from a
    near one, so a failure from far away says little about pull-in.

!!! note "Rounding in the forces"
    Forces of the field method carry rounding of its dense solutions, which
    grows with the size of a cluster. A structure with a very soft motion
    turns that into steps above a tight $\tau_s$, which is what the second way
    of stopping is for.

!!! note "Cost of an iterate"
    Every iterate costs one evaluation of the capacitance model at
    displaced bodies. That is a fraction of a second with plates alone or
    with the plates method, and with the field method as long as the model
    at rest took for the clusters whose bodies have moved in the plane; a
    cluster with none moved keeps the field solved at rest, since rise and
    tilt do not enter the sidewalls.

!!! note "Limits that carry over"
    Beams stay linear at any displacement, and the limits of the
    capacitance model at displaced bodies carry over: tilt of a plate to
    second order, and in-plane travel of a plate to first.

!!! note "Loose bodies"
    A body that no chain of beams joins to the anchor has no stiffness, and
    nothing balances a force on it. Such a body, for instance an island
    left on a neck too thin to qualify as a beam, is kept at rest and
    listed with the result. (120) to (124) are taken over the other bodies.

## Symbols on this page

| Symbol | Meaning | In the code |
| --- | --- | --- |
| $\mathbf{v}^\text{net}$, $\mathbf{Q}_f$ | potentials of driven nets and charges of floating ones | `potentials` and `charges` of `Structure.equilibrium`; `Equilibrium.floating` |
| $\mathbf{f}_0$ | constant mechanical loads | `loads` of `Structure.equilibrium` |
| $\mathbf{q}$ | displacement at equilibrium | `Equilibrium.displacement`, `Equilibrium.dofs` |
| $\mathbf{q}_0$ | start of the iteration | `start` of `Structure.equilibrium` |
| $\tau_s$ | tolerance on the last step | `StaticOptions.tolerance` |
| $\mathbf{K}_t$ | tangent stiffness, with the part a floating net adds | `Equilibrium.stiffness` |
| $\lambda_\text{min}$ | squared angular frequency of the softest motion under bias | `Equilibrium.softest`, `Equilibrium.stable` |
| $\mu_\text{min}$ | share of beam stiffness left where the field softens most | `Equilibrium.margin` |
| $\chi$ | clearance: narrowest gap against its width at rest | `CapacitanceModel.clearance`, `CapacitanceModel.closest` |
| $V^*$ | pull-in voltage | found by the caller's search, as the usage page shows |
