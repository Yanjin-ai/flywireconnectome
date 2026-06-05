# A Conserved Sensorimotor Circuit Identified Across Three Drosophila Connectomes

**Dataset:** BANC (brain & nerve cord, ♀) × FAFB (adult brain, ♀) × MANC (nerve cord, ♂)  
**Circuit size:** N = 75 neurons · 6 verified conserved synaptic edges  
**Method:** Neuron-level MCIS with NBLAST-matched cross-connectome correspondence

---

## Circuit Visualization

| Figure | Description |
|--------|-------------|
| ![Circuit layouts](figures/figure1_circuit_layouts.png) | **Fig. 1** — Three force-directed layouts of the 75-neuron circuit. Gold edges: the 6 synaptic connections verified identical across all three connectomes. Node colour encodes neuron class (red = descending, blue = ascending). |
| ![Composition](figures/figure2_composition.png) | **Fig. 2** — Neuron class composition, neurotransmitter profile, motor target regions, degree distribution, and cross-dataset edge count comparison. |
| ![Hub circuit](figures/figure3_hub_circuit.png) | **Fig. 3** — Hub neurons (large nodes) connected by the 6 conserved edges. Labels show cell-type identifiers (e.g. DNp63, DNge059). |
| ![Dimorphism](figures/figure4_dimorphism_nt.png) | **Fig. 4** — Sexual dimorphism status and neurotransmitter × neuron-class breakdown. 91% of neurons are sexually isomorphic. |

---

## Biological Significance

### What is this circuit?

The 75-neuron circuit is dominated by **Descending Neurons** (DN, n = 46) and **Ascending Neurons** (AN, n = 22), with a minority of sensory-ascending neurons (n = 4) and sensory-descending neurons (n = 3). These neuron classes form the anatomical and functional **sensorimotor interface** of the insect central nervous system — the bidirectional information highway connecting the decision-making brain to the movement-executing ventral nerve cord (VNC).

Descending neurons carry motor commands from the brain to the VNC; ascending neurons relay proprioceptive and mechanosensory feedback from the body back to the brain. Their tight co-conservation in a single isomorphic subgraph — across a female brain (FAFB), a female brain-plus-cord (BANC), and a male nerve cord (MANC) — is a direct signature of the **sensorimotor bottleneck** postulated by Pospisil et al. (2024).

### Neurotransmitter profile

| Neurotransmitter | Count | Interpretation |
|-----------------|-------|----------------|
| Acetylcholine   | 47 (63%) | Excitatory fast transmission; dominant in insect motor control |
| GABA            | 15 (20%) | Inhibitory; consistent with gain-control roles in DN populations |
| Glutamate       |  9 (12%) | Mixed excitatory/inhibitory depending on receptor type |
| Serotonin       |  3  (4%) | Neuromodulatory; state-dependent modulation of locomotor circuits |
| Dopamine        |  1  (1%) | Reward/motivational modulation |

### Motor target regions

The cns_network labels of constituent neurons reveal that this circuit targets primarily:
- **Leg VNC** (n = 25): locomotion control
- **Dorsal VNC** (n = 15): wing and flight motor neurons
- **Flange median bundle** (n = 8): descending tract mediating whole-body coordination
- **Abdominal VNC** (n = 7): abdominal and reproductive motor control

Collectively, this points to a **multi-effector motor coordination hub** that integrates commands for walking, flight, and postural control simultaneously.

### Sexual conservation

**91% of neurons (68/75) are annotated as sexually isomorphic** (identical in male and female), with only 7 neurons showing sex-specific differences. This is highly significant: the core wiring of the circuit is preserved across sexes despite substantial sexual dimorphism elsewhere in the fly brain (Berg et al. 2025 report ~4.8% dimorphism in central brain cell types). The isomorphic majority suggests that the circuit encodes computations essential to both sexes — basic locomotor coordination is sex-independent.

### The 6 conserved edges

The 6 edges verified identical across BANC, FAFB, and MANC represent connections that survive all sources of biological and technical variability:
- Individual-to-individual variation (different animals)
- Sex differences (♀ FAFB/BANC vs. ♂ MANC)
- Dataset reconstruction differences (EM segmentation artefacts)

These edges therefore reflect **genetically encoded, functionally obligate synaptic connections** — the hardwired backbone of sensorimotor signal flow.

---

## Hypothesis

> **The identified 75-neuron circuit constitutes the genetically canalized sensorimotor bottleneck of the Drosophila CNS.** Its inter-sex, inter-individual invariance reflects a developmental constraint on the brain-body interface: descending motor commands and ascending sensory feedback must flow through a conserved scaffold, while higher-order circuits (mushroom body, central complex) remain plastic and sex-dimorphic.

This predicts:
1. Silencing any hub DN in this circuit (e.g. DNp63, DNpe016) should impair multiple motor programs simultaneously — not just one behaviour.
2. The 6 conserved edges should have higher synapse counts (stronger connections) than average, making them robust to stochastic noise.
3. Orthologous circuits should be identifiable in other insects (e.g. *Manduca sexta*, *Apis mellifera*) wherever connectome data become available.

---

## Key Literature

1. **Dorkenwald et al. (2024)** — Neuronal wiring diagram of an adult brain. *Nature* 634, 123–138. [doi:10.1038/s41586-024-07558-y](https://doi.org/10.1038/s41586-024-07558-y)
2. **Schlegel et al. (2024)** — Whole-brain annotation and multi-connectome cell typing of *Drosophila*. *Nature* 634, 139–152. [doi:10.1038/s41586-024-07686-5](https://doi.org/10.1038/s41586-024-07686-5) — *Source of NBLAST cross-connectome matching methodology.*
3. **Pospisil et al. (2024)** — The fly connectome reveals a path to the effectome. *Nature* 634, 234–242. [doi:10.1038/s41586-024-07982-0](https://doi.org/10.1038/s41586-024-07982-0) — *Sensorimotor bottleneck concept.*
4. **Berg et al. (2025)** — Sexual dimorphism in the complete connectome of the *Drosophila* male central nervous system. *bioRxiv*. [doi:10.1101/2025.10.09.680999](https://doi.org/10.1101/2025.10.09.680999) — *Cross-sex conservation baseline.*
5. **Bates et al. (2025)** — Distributed control circuits across a brain-and-cord connectome. *bioRxiv*. [doi:10.1101/2025.07.31.667571](https://doi.org/10.1101/2025.07.31.667571) — *BANC-centred comparative analysis; motivates choice of triplet.*
6. **Takemura et al. (2024)** — A connectome of the male *Drosophila* ventral nerve cord. *eLife* 13, e97769. [doi:10.7554/eLife.97769](https://doi.org/10.7554/eLife.97769) — *MANC dataset.*
7. **Witvliet et al. (2021)** — Connectomes across development reveal principles of brain maturation. *Nature* 596, 257–261. [doi:10.1038/s41586-021-03778-8](https://doi.org/10.1038/s41586-021-03778-8) — *Methodology for connectome stereotypy analysis.*
8. **Milo et al. (2002)** — Network motifs: simple building blocks of complex networks. *Science* 298, 824–827. [doi:10.1126/science.298.5594.824](https://doi.org/10.1126/science.298.5594.824) — *Circuit motif framework.*
