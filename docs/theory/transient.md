# Transient

Under drive that is large, switched, or fed through resistance, the structure
is followed in time: the motion of its bodies, the charge on the nodes that
resistance sets apart, and the potentials that follow from both. This page
gives the equations of motion and of the network, the integrator and its
tolerances, the current of each source, the energy ledger that checks the
integration, the start of a run, and the contact event that ends one early.
`Structure.transient` runs it and returns a `Transient`.

Under large or switched drive, through resistance, the structure is followed
in time. The unknowns are the displacement $\mathbf{q}$ and velocity
$\mathbf{u}$ of the bodies that beams join to the anchor, and the charge
$\mathbf{Q}_d$ on each node that a source reaches only through resistance.
With the force of (120) at the potentials of the moment and constant
mechanical loads $\mathbf{f}_0$,

$$
\dot{\mathbf{q}} = \mathbf{u}, \qquad
\mathbf{M}^\text{full} \dot{\mathbf{u}} + \mathbf{D}\, \mathbf{u} + \mathbf{K}^\text{full} \mathbf{q} = \mathbf{f}(\mathbf{q}, \mathbf{v}) + \mathbf{f}_0, \qquad
\dot{\mathbf{Q}}_d = \mathbf{i}_d. \tag{139}
$$

Damping is that of one quality factor $Q$ for every mode of the structure
without potentials, the modes $\boldsymbol{\Phi}$ and frequencies
$\omega_n$ of (38):

$$
\mathbf{D} = \mathbf{M}^\text{full} \boldsymbol{\Phi}\, \operatorname{diag}\Bigl( \frac{\omega_n}{Q} \Bigr)\, \boldsymbol{\Phi}^\top \mathbf{M}^\text{full}. \tag{140}
$$

A fixed matrix keeps the ledger of (145) exact. Under bias the modes
soften: a mode that keeps its shape has the quality factor $Q\, \omega'_n
/ \omega_n$ at its frequency $\omega'_n$, slightly below $Q$, and where
bias mixes the shapes $\mathbf{D}$ couples the biased modes, so the damping
of (132) about an equilibrium and that of a transient agree only without
bias.

### Network

Nodes are those of (130), with each voltage source in place of a port: it
drives the node of its electrode. The net of an electrode that a current
source or an element of the circuit around the structure reaches joins the
nodes as the net of a port does. An ideal source, without resistance, holds
its node $h$ at its waveform $e_h(t)$. A source of resistance $R_k$ feeds its
node $k$ through it, and its node is one of the nodes $d$ that carry charge
as a state, with every other node of a driven net. A floating net is one
node $f$ that keeps its charge $\mathbf{Q}_f$, and a net that nothing reaches
is held at zero volts on every electrode: it is ground. Resistors of the
circuit add to $\mathbf{G}$ of (130) as links do, a resistor to ground on the
diagonal alone. A switch $s(t)$, one while closed and zero while open,
multiplies the branch it is on, a source's or a resistor's, so
$\mathbf{G}(t)$ and the feeds change at its times. The current into a node
$d$ is

$$
i_d = -\sum_l G_{dl}(t)\, v_l + s_d(t)\, \frac{e_d - v_d}{R_d} + j_d(t), \tag{141}
$$

with the second term only at a node that a source with resistance feeds and
$j_d$ the current that current sources drive into it. With $U$ the nodes $d$
and $f$ together, $\mathbf{C}$ is the node capacitance of (130) at
$\mathbf{q}$ with the capacitors $C_c$ of the circuit added, and the charges
fix the potentials:

$$
\mathbf{C} = \mathbf{T}^\top \mathbf{C}^\text{el}(\mathbf{q})\, \mathbf{T} + \sum_c C_c\, (\mathbf{e}_k - \mathbf{e}_l)(\mathbf{e}_k - \mathbf{e}_l)^\top, \qquad
\mathbf{C}_{UU}\, \mathbf{v}_U = \begin{bmatrix} \mathbf{Q}_d \\ \mathbf{Q}_f \end{bmatrix} - \mathbf{C}_{Uh}\, \mathbf{e}_h(t), \tag{142}
$$

where $\mathbf{C}^\text{el}$ is the matrix between electrodes, capacitor $c$
lies between nodes $k$ and $l$, and $\mathbf{e}_l$ is zero for a capacitor
to ground.

Charge rather than potential is the state because charge is what
resistance limits: a source that changes moves the potentials of other
nodes at once through capacitance, but no charge through a resistor.
Solving (142) at every evaluation leaves (139) an ordinary differential
system without constraints.

### Integration

Resistance and capacitance of beams and sources give time constants that
are often far shorter than mechanical periods, so the system is stiff. It
is integrated by the fifth-order Radau IIA method, which is L-stable,
with error control and the Jacobian of (139) and (142) in closed form; the
rows of the three energy states of (145) are left zero, since nothing feeds
back from them, which slows the integrator's own Newton iteration at most.
From (129) at node
level, $\boldsymbol{\Gamma}$ with rows over nodes, and with $\tilde{\mathbf{G}}$
the conductance $\mathbf{G}(t)$ of (141) with $s_k / R_k$ added on the
diagonal at nodes fed through resistance and $\mathbf{E}_d$ the columns of
the nodes $d$,

$$
\frac{\partial \mathbf{f}}{\partial \mathbf{q}} = -\mathbf{K}_e - \boldsymbol{\Gamma}_U^\top \mathbf{C}_{UU}^{-1} \boldsymbol{\Gamma}_U, \quad
\frac{\partial \mathbf{f}}{\partial \mathbf{Q}_d} = \boldsymbol{\Gamma}_U^\top \mathbf{C}_{UU}^{-1} \mathbf{E}_d, \quad
\frac{\partial \mathbf{i}_d}{\partial \mathbf{q}} = \tilde{\mathbf{G}}_{dU}\, \mathbf{C}_{UU}^{-1} \boldsymbol{\Gamma}_U, \quad
\frac{\partial \mathbf{i}_d}{\partial \mathbf{Q}_d} = -\tilde{\mathbf{G}}_{dU}\, \mathbf{C}_{UU}^{-1} \mathbf{E}_d, \tag{143}
$$

with $\mathbf{K}_e$ of (116). The first is (128) extended to every node
whose charge, not potential, is held over the instant. Current sources and
capacitors of the circuit add nothing to it beyond $\tilde{\mathbf{G}}$ and
$\mathbf{C}$. Integration restarts at the breakpoints of each waveform, where
its slope changes, and at the times of each switch.

### Current of a source

A source with resistance drives the current $s_k (e_k - v_k) / R_k$. An
ideal source drives what its node takes, the change of its charge and the
current into its links and resistors, less what current sources drive into
it:

$$
i_h = \bigl( \dot{\mathbf{C}} \mathbf{v} + \mathbf{C} \dot{\mathbf{v}} + \mathbf{G} \mathbf{v} \bigr)_h - j_h, \qquad
\mathbf{C}_{UU}\, \dot{\mathbf{v}}_U = \begin{bmatrix} \dot{\mathbf{Q}}_d \\ \mathbf{0} \end{bmatrix} - \bigl( \dot{\mathbf{C}} \mathbf{v} \bigr)_U - \mathbf{C}_{Uh}\, \dot{\mathbf{e}}_h, \qquad
\dot{\mathbf{C}} = \sum_k \frac{\partial \mathbf{C}}{\partial q_k}\, u_k, \tag{144}
$$

which needs the slope $\dot{e}_h$ of its waveform.

### Energy

The energy stored in the bodies, the beams, the field, and the capacitors
of the circuit, and the energy taken and supplied,

$$
\begin{aligned}
W &= \tfrac{1}{2} \mathbf{u}^\top \mathbf{M}^\text{full} \mathbf{u} + \tfrac{1}{2} \mathbf{q}^\top \mathbf{K}^\text{full} \mathbf{q} + \tfrac{1}{2} \mathbf{v}^\top \mathbf{C}\, \mathbf{v}, \\
W_D &= \int \mathbf{u}^\top \mathbf{D}\, \mathbf{u}\, dt, \qquad
W_R = \int \Bigl( \mathbf{v}^\top \mathbf{G}(t)\, \mathbf{v} + \sum_k s_k\, \frac{(e_k - v_k)^2}{R_k} \Bigr) dt, \qquad
W_S = \int \Bigl( \sum_s e_s\, i_s + \sum_d v_d\, j_d + \mathbf{f}_0^\top \mathbf{u} \Bigr) dt,
\end{aligned} \tag{145}
$$

balance exactly, $W(t) - W(t_0) + W_D + W_R = W_S$, with the sum over every
voltage source $s$ in $W_S$, $\mathbf{C}$ of (142), and the potentials of
every node in $W$. The three
integrals are integrated as states beside (139), so what remains of the
balance is the error of integration.

### Tolerances and start

Each step keeps its error below $\tau_t$ times each state plus an absolute
part, in units that make states of different kinds comparable:

$$
\tau_t\, g_0, \quad \tau_t\, \frac{g_0}{\sqrt{A_B}}, \quad \tau_t\, g_0\, \omega_1, \quad \tau_t\, \frac{g_0}{\sqrt{A_B}}\, \omega_1, \quad \tau_t\, C_{dd} V_0, \quad \tau_t\, \tfrac{1}{2} V_0^2 \sum_k C_{kk}, \tag{146}
$$

for translations, rotations of a body of area $A_B$, its device silicon or
its handle silicon when it has none, with $\sqrt{A_B}$ taken no smaller
than $g_0$, translational and angular velocities, charges, and energies,
with $g_0$ the gap of the buried oxide, $\omega_1$ the lowest
frequency of (38), and $V_0$ the largest source value at the output times
and breakpoints, at least one volt. The run starts at $\mathbf{q}_0$, rest
or the displacement given, at zero velocity, with the potentials that
carry no current at $t_0$. Nodes $D$ that resistance ties, closed at $t_0$,
to a source or to ground take

$$
\tilde{\mathbf{G}}_{DD}\, \mathbf{v}_D = \Bigl( \frac{s_d\, e_d(t_0)}{R_d} + j_d(t_0) \Bigr) - \tilde{\mathbf{G}}_{Dh}\, \mathbf{e}_h(t_0), \qquad
\sum_{d \in g} \bigl( \mathbf{C} \mathbf{v} \bigr)_d = 0, \qquad
\bigl( \mathbf{C} \mathbf{v} \bigr)_f = Q_f, \qquad
\mathbf{Q}_d = \bigl( \mathbf{C} \mathbf{v} \bigr)_d, \tag{147}
$$

where the first term on the right is zero at nodes without a source. Nodes
that resistance joins to neither form groups $g$ that resistance joins
among themselves, each at one potential and with no charge in all: a pad
behind a capacitor alone, or behind a switch open at $t_0$. Floating nets
keep their charge, and the last two parts are solved together.

### Contact

The run ends when the clearance $\chi$ of (126) falls to $\chi_c$. The
moment is found on the continuous extension of the last step, and the
result reports it last, with the electrodes across the closing gap.

!!! note "Switches, current sources, and what is not built"
    A switch changes at its times, and an output at one of them shows the
    branch changed. A switched source needs resistance: an ideal source
    holds its node, and opening it would turn a held potential into a
    charge in the middle of a run. A current source drives its node from
    ground. Inductors and slew limits are not built: electrical time
    constants are far shorter than mechanical ones, and a driver's slew is a
    ramp in its waveform.

!!! note "Waveforms must be continuous"
    A waveform must be continuous: a step is a ramp of the rise time it has.
    A drive that switches in no time would deliver its energy in no time
    through an ideal source, which (144) cannot follow, and no real driver
    does so. A waveform that is a `PiecewiseLinear` gives its breakpoints
    and slopes; another callable is differentiated by central differences
    for (144) and should be continuous between breakpoints it lists.

!!! note "Accuracy of outputs"
    Outputs come from the continuous extension of the step that spans them,
    which is less accurate than the ends of the steps. Every output
    evaluates the capacitance model once more.

!!! note "Every degree of freedom integrated"
    Every degree of freedom of the bodies that beams join to the anchor is
    integrated; no mode is left out. Bodies that no beam holds stay at
    rest. A floating net is one node, as in small signal.

!!! note "Cost of a step"
    Each step evaluates the capacitance model at the stages of the method
    and its Newton iterations, about ten times. That is fractions of a
    millisecond with plates alone, a few milliseconds with the plates
    method, and as long as the solve at rest with the field method, which
    makes long runs on the field method slow.

!!! note "Steps into contact"
    A trial step that would bring two bodies into contact is retried with
    a step four times shorter, up to eight times; the run then ends at the
    last accepted step with a contact event.

## Symbols on this page

| Symbol | Meaning | In the code |
| --- | --- | --- |
| $\mathbf{q}$, $\mathbf{u}$ | displacement and velocity of the bodies | `Transient.displacement`, `Transient.velocity` |
| $\mathbf{Q}_d$, $\mathbf{v}$ | charges on nodes behind resistance and potentials of every electrode | `Transient.charge`, `Transient.potentials` |
| $\mathbf{D}$, $Q$ | damping matrix and the quality factor that sets it | `TransientOptions.quality_factor` |
| $e_k(t)$, $R_k$, $s_k(t)$ | waveform, resistance, and switch of a voltage source | `Source.waveform`, `Source.resistance`, `Source.switch` |
| $j_d(t)$ | current of a current source | `CurrentSource.waveform` |
| $\mathbf{G}(t)$, $C_c$ | conductance of links and resistors, capacitors of the circuit | `Resistor`, `Capacitor` in `circuit` |
| $i_h$, $i_k$ | current of each source | `Transient.current` |
| $W$, $W_D$, $W_R$, $W_S$ | stored, damped, resistive, and supplied energy | `Transient.energy` and its `imbalance` |
| $\tau_t$ | tolerance of the integrator | `TransientOptions.tolerance`, `TransientOptions.max_step` |
| $\chi_c$ | clearance at which the run ends | `TransientOptions.contact_clearance`, `Transient.contact`, `Transient.clearance` |
