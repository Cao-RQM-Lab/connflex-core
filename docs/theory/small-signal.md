# Small signal

About a stable equilibrium, small signals at a few ports move the bodies by
small amounts and draw small currents, and the structure is linear in both.
This page derives that linear model: the transduction between potentials and
motion, the network of nodes that resistance in the links sets apart, the
modes under bias with their damping, the admittance between ports, the
motional branches of each mode, the range in which the linearization holds,
and the SPICE netlist that represents it. `Structure.small_signal` builds it
from an `Equilibrium` and returns a `SmallSignal`.

About an equilibrium $\mathbf{q}_0$ at electrode potentials $\mathbf{v}_0$,
small changes $\delta \mathbf{v}(t)$ of the potentials move the bodies by
small $\delta \mathbf{q}(t)$. To first order in both, the force of (120) and
the charge $\mathbf{Q} = \mathbf{C} \mathbf{v}$ change by

$$
\delta \mathbf{f} = -\mathbf{K}_e\, \delta \mathbf{q} + \boldsymbol{\Gamma}^\top \delta \mathbf{v}, \qquad
\delta \mathbf{Q} = \mathbf{C}\, \delta \mathbf{v} + \boldsymbol{\Gamma}\, \delta \mathbf{q}, \qquad
\Gamma_{ak} = \sum_b \frac{\partial C_{ab}}{\partial q_k}\, v_{0,b}, \tag{129}
$$

with $\mathbf{C}$, its gradient, and $\mathbf{K}_e$ of (116) taken at
$\mathbf{q}_0$. $\boldsymbol{\Gamma}$ has one row per electrode and one
column per degree of freedom; summed over the electrodes of each net it is
the matrix of (128). The same $\boldsymbol{\Gamma}$ turns potentials into
force and motion into charge, so a port senses a motion with the
coefficient by which it drives it, and the model is reciprocal.

### Network

Signals carry current, so links carry voltage and the electrodes of a net
no longer share a potential. A source holds each port, an electrode named
by the caller, at its potential in the equilibrium plus the signal. A net
with a potential and no port is held at every electrode and drops out. A
floating net is one node. Electrodes that links without resistance join
are one node as well, and the nodes left are the ports and the electrodes
of their nets that links with resistance set apart from them: the signal
nodes $s$, ports $p$ first and inner nodes $i$ after them. Floating nodes
are $f$. With $T_{ak} = 1$ when electrode $a$ belongs to node $k$ and zero
otherwise,

$$
\mathbf{C}^\text{node} = \mathbf{T}^\top \mathbf{C}\, \mathbf{T}, \qquad
\boldsymbol{\Gamma}^\text{node} = \mathbf{T}^\top \boldsymbol{\Gamma}, \qquad
\mathbf{G} = \sum_\text{links} \frac{1}{R}\, (\mathbf{e}_k - \mathbf{e}_l)(\mathbf{e}_k - \mathbf{e}_l)^\top, \tag{130}
$$

where the sum runs over links of resistance $R$ between signal nodes $k
\ne l$, with $R$ of (41) for beams.

A floating node keeps its charge, $\delta Q_f = 0$, which by (129) sets its
potential. Eliminating it,

$$
\delta \mathbf{v}_f = -\mathbf{C}_{ff}^{-1} \bigl( \mathbf{C}_{fs}\, \delta \mathbf{v}_s + \boldsymbol{\Gamma}_f\, \delta \mathbf{q} \bigr), \qquad
\mathbf{C}' = \mathbf{C}_{ss} - \mathbf{C}_{sf}\, \mathbf{C}_{ff}^{-1} \mathbf{C}_{fs}, \qquad
\boldsymbol{\Gamma}' = \boldsymbol{\Gamma}_s - \mathbf{C}_{sf}\, \mathbf{C}_{ff}^{-1} \boldsymbol{\Gamma}_f, \tag{131}
$$

in blocks of (130). The charge on the signal nodes is then $\delta
\mathbf{Q}_s = \mathbf{C}' \delta \mathbf{v}_s + \boldsymbol{\Gamma}' \delta
\mathbf{q}$, and the force $\boldsymbol{\Gamma}'^\top \delta \mathbf{v}_s$
less $\boldsymbol{\Gamma}_f^\top \mathbf{C}_{ff}^{-1} \boldsymbol{\Gamma}_f\,
\delta \mathbf{q}$, which is the part that (128) adds to the tangent
stiffness. The stiffness of the linear model is therefore
$\mathbf{K}_t$, that of the equilibrium.

### Modes and damping

The motion is taken in the modes of the structure under bias,

$$
\mathbf{K}_t\, \boldsymbol{\phi}_n = \omega_n^2\, \mathbf{M}^\text{full} \boldsymbol{\phi}_n, \quad
\boldsymbol{\phi}_n^\top \mathbf{M}^\text{full} \boldsymbol{\phi}_m = \delta_{nm}, \qquad
\ddot{a}_n + \frac{\omega_n}{Q}\, \dot{a}_n + \omega_n^2\, a_n = \boldsymbol{\eta}_n^\top \delta \mathbf{v}_s, \qquad
\boldsymbol{\eta}_n = \boldsymbol{\Gamma}' \boldsymbol{\phi}_n, \tag{132}
$$

with $\delta \mathbf{q} = \sum_n \boldsymbol{\phi}_n a_n$, $\mathbf{M}^\text{full}$
of (34), and one quality factor $Q$ for every mode: the damping matrix
$\mathbf{M}^\text{full} \boldsymbol{\Phi}\, \operatorname{diag}(\omega_n / Q)\,
\boldsymbol{\Phi}^\top \mathbf{M}^\text{full}$, which no mode couples to
another. $\boldsymbol{\eta}_n$ is the generalized force on mode $n$ per
volt on each signal node. $\mathbf{K}_t$ is softer than
$\mathbf{K}^\text{full}$ where the field pulls toward a closing gap, so
$\omega_n$ lies below the frequency of (38) without bias. For a drive
$\delta \mathbf{v}_s\, e^{j \omega t}$,

$$
\delta \mathbf{q} = \sum_n \boldsymbol{\phi}_n\, \frac{\boldsymbol{\eta}_n^\top \delta \mathbf{v}_s}{\Delta_n(\omega)}, \qquad
\Delta_n(\omega) = \omega_n^2 - \omega^2 + j \omega\, \frac{\omega_n}{Q}. \tag{133}
$$

### Admittance

The current that enters a signal node from outside charges it and flows
on through its links:

$$
\delta \mathbf{i}_s = \mathbf{A}(\omega)\, \delta \mathbf{v}_s, \qquad
\mathbf{A}(\omega) = \mathbf{G} + j \omega \Bigl( \mathbf{C}' + \sum_n \frac{\boldsymbol{\eta}_n \boldsymbol{\eta}_n^\top}{\Delta_n(\omega)} \Bigr). \tag{134}
$$

No current enters an inner node from outside, so

$$
\delta \mathbf{v}_i = -\mathbf{A}_{ii}^{-1} \mathbf{A}_{ip}\, \delta \mathbf{v}_p, \qquad
\mathbf{Y}(\omega) = \mathbf{A}_{pp} - \mathbf{A}_{pi}\, \mathbf{A}_{ii}^{-1} \mathbf{A}_{ip}, \tag{135}
$$

and $\mathbf{Y}$ is the admittance between ports: entry $(p, r)$ is the
current into port $p$ per volt at port $r$ with the other ports held.
It is symmetric, by reciprocity, and its Hermitian part is positive
semidefinite, since resistors and damping only take energy.

### Motional branches

With every other node held, mode $n$ adds $j \omega\, \eta_{ns}^2 /
\Delta_n$ to the admittance at node $s$, which is that of a series branch
to ground:

$$
\frac{j \omega\, \eta_{ns}^2}{\Delta_n(\omega)} = \Bigl( j \omega L_{ns} + R_{ns} + \frac{1}{j \omega C_{ns}} \Bigr)^{-1}, \qquad
L_{ns} = \frac{1}{\eta_{ns}^2}, \quad
C_{ns} = \frac{\eta_{ns}^2}{\omega_n^2}, \quad
R_{ns} = \frac{\omega_n}{Q\, \eta_{ns}^2}. \tag{136}
$$

A port that is the only node, near a resonance that couples to it alone,
is $\mathbf{C}'$ in parallel with this branch, the Butterworth-Van Dyke
circuit. For a plate of mass $m$ on a spring $k_t$ under bias $V$ across a
gap $g$ of area $A$, $\eta^2 = \Gamma^2 / m$ with $\Gamma = \varepsilon A V
/ g^2$, so $L = m / \Gamma^2$, $C = \Gamma^2 / k_t$, and $R = \sqrt{k_t
m} / (Q\, \Gamma^2)$.

### Range of the linearization

Two expansions in (129) fail first. The force between electrodes $a$ and
$b$ goes as $(v_a - v_b)^2$, and for a swing $\delta$ of that difference
about a bias $\Delta v$ the first order drops $\delta^2$ against $2 \Delta
v\, \delta$: a relative error of half the ratio $\delta / \Delta v$, which
appears at twice the drive frequency. Capacitance across a gap $g$ goes as
$1 / g$ and its gradient as $1 / g^2$, so a motion that closes the
narrowest gap by a fraction $r$ changes the coefficients of (129) by about
$2 r$. For a drive at one frequency the two ratios are

$$
r_v = \max_{(a, b)} \frac{|\delta v_a - \delta v_b|}{|v_{0,a} - v_{0,b}|}, \qquad
r_g = 1 - \frac{\min_\pm \chi(\mathbf{q}_0 \pm \mathbf{p})}{\chi(\mathbf{q}_0)}, \tag{137}
$$

the first over pairs of electrodes whose capacitance changes with motion,
with potentials of every electrode from (135) and (131), and infinite for a
pair without bias that swings. In the second, $\chi$ is the clearance of
(126), and $\mathbf{p}$ is the displacement at the moment of the cycle
$\operatorname{Re}\bigl( \delta \mathbf{q}\, e^{j \theta} \bigr)$ that
stores most strain energy $\mathbf{p}^\top \mathbf{K}^\text{full}
\mathbf{p}$. Either above $\tau_l$ draws a warning; both are reported, and
what error a design tolerates is the caller's judgment.

### Netlist

A mode that couples to several nodes is one series loop that all of them
drive and all of them feel. With $\sigma_n = \max_s |\eta_{ns}|$, a
loop of $L_n = 1 / \sigma_n^2$, $R_n = \omega_n / (Q \sigma_n^2)$, and $C_n
= \sigma_n^2 / \omega_n^2$ in series with sources of gains $\eta_{ns} /
\sigma_n$ carries

$$
I_n = \frac{\sum_s (\eta_{ns} / \sigma_n)\, \delta v_s}{j \omega L_n + R_n + 1 / (j \omega C_n)}
= \frac{j \omega\, \sigma_n\, \boldsymbol{\eta}_n^\top \delta \mathbf{v}_s}{\Delta_n(\omega)}, \qquad
\delta i_s = \frac{\eta_{ns}}{\sigma_n}\, I_n, \tag{138}
$$

which is the term of mode $n$ in (134). In SPICE the sources in the loop
are voltage-controlled voltage sources in series, each controlled by one
node, and the current returns to each node through a current-controlled
current source. Capacitors join signal nodes $k$ and $l$ with
$-C'_{kl}$ and each node to ground with its row sum $\sum_l C'_{kl}$,
which is its capacitance to held nets, and resistors join them with
$-1 / G_{kl}$. Ground is the potential of nets held without a port. With
one node, $\sigma_n = |\eta_{ns}|$ and the loop is (136) with gains of
$\pm 1$.

!!! note "Ports are electrodes"
    Ports are electrodes rather than nets because a source drives a net at
    a pad, and inner electrodes behind beams of resistance see it later.
    A net held without a port is held at every electrode, which neglects
    its own resistance; making it a port restores it. Ports are given as
    electrode objects, never as indices: nets and electrodes are both
    numbered from zero, and an index of one would pass for the other.

!!! note "Floating nets"
    A floating net is one node, as it is in equilibrium. Charge that moves
    between its electrodes through resistance is left out.

!!! note "Nodes joined without resistance"
    Links without resistance join electrodes into one node: bridges, ideal
    vias, ideal wires, and beams when the device material has no
    resistivity. A net without resistivity is therefore one node, and two
    ports on it are an error, since two sources would hold one node.

!!! note "A stable equilibrium is required"
    The equilibrium must be stable. An unstable one has a mode with
    $\omega_n^2 < 0$, which no branch represents.

!!! note "Modes and the netlist subset"
    Modes come from a dense eigensolution over the bodies that beams join
    to the anchor, all of them; loose bodies stay at rest. A netlist may
    take a subset, and a mode left out takes with it its compliance below
    resonance, $\eta_{ns}^2 / \omega_n^2$ of capacitance at node $s$.

!!! note "Couplings set to zero"
    A coupling smaller than $D \epsilon\, \bigl( \sum_n \eta_{ns}^2
    \bigr)^{1/2}$, with $D$ the degrees of freedom kept and $\epsilon$ the
    machine epsilon, is set to zero. That sum is $\boldsymbol{\Gamma}'_s
    \mathbf{M}^{\text{full}\,-1} \boldsymbol{\Gamma}_s'^\top$ for any set of
    mass-normalized modes, and couplings below it by this factor are
    rounding of the mode shapes: modes that symmetry decouples from a node
    get no branch there.

!!! note "Damping"
    Damping enters only through $Q$. Squeeze film, gas, anchor, and
    thermoelastic loss are not computed.

!!! note "Zero frequency and the gap ratio"
    Frequencies may be zero, where $\mathbf{A} = \mathbf{G}$ and
    $\mathbf{Y}$ is the conductance between ports of one net. The gap
    ratio of (137) evaluates the capacitance model at both extremes of the
    motion, which costs two evaluations of it.

## Symbols on this page

| Symbol | Meaning | In the code |
| --- | --- | --- |
| $\boldsymbol{\Gamma}$ | transduction: charge per motion and force per volt | `CapacitanceModel.charge_gradient` at the equilibrium |
| $p$, $i$, $f$ | ports, inner nodes, and floating nodes | `ports` of `Structure.small_signal`; `SmallSignal.nodes`, ports first |
| $\mathbf{C}'$, $\mathbf{G}$ | node capacitance after eliminating floating nodes, and node conductance | `SmallSignal.capacitance`, `SmallSignal.conductance` |
| $\omega_n$, $\boldsymbol{\phi}_n$ | frequency and shape of a mode under bias | `SmallSignal.frequencies_hz`, `SmallSignal.modes` |
| $Q$ | quality factor of every mode | `SmallSignalOptions.quality_factor` |
| $\boldsymbol{\eta}_n$ | force on mode $n$ per volt on each node | `SmallSignal.transduction` |
| $\mathbf{Y}(\omega)$ | admittance between ports | `SmallSignal.admittance` |
| $\delta \mathbf{q}$ | complex amplitude of the motion | `SmallSignal.displacement` |
| $L_{ns}$, $C_{ns}$, $R_{ns}$ | motional branch of mode $n$ at node $s$ | `SmallSignal.motional` |
| $r_v$, $r_g$, $\tau_l$ | voltage and gap ratios of a drive and the tolerance that warns | `SmallSignal.linearity`, `SmallSignalOptions.linear_tolerance` |
| $\sigma_n$ | scale of the loop of mode $n$ in the netlist | `SmallSignal.netlist` |
