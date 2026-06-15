# Results and Discussion Paragraphs for the Recommended Figure Sequence

## Writing Diagnosis

Current stage: Stage 3, Mechanism Construction, with sentence-level claim calibration.

Core contradiction: The manuscript should not read as "RT outputs were connected to an SE model." The stronger story is that dissolution changes hydraulic connectivity and pore-fluid conductivity, and their competition controls the interfacial electrokinetic current imbalance that is transferred into $L(\omega)$, $\sigma(\omega)$, $R_E/T_{\mathrm{TM}}$, and finite-offset waveform behavior.

Currently forbidden: claiming experimental validation, claiming field detectability, or overstating monitoring capability beyond the current simulations.

## Figure 1

Figure 1 defines the two time axes that organize the analysis. Reactive transport evolves over dissolution time, during which porosity, permeability, tortuosity, $\mathrm{H}^+$ concentration, and pore-fluid conductivity change gradually. For each dissolution state, these variables are mapped into dynamic electrokinetic properties and then into a waveform calculated over the much shorter acoustic response time. This separation prevents the waveform time series from being interpreted as chemical evolution, and it clarifies the mechanism chain used throughout the study: dissolution modifies hydraulic and electrical properties, these properties reshape $L(\omega)$ and $\sigma(\omega)$, and the resulting interfacial current imbalance controls the finite-offset seismoelectric response.

## Figure 2

Figure 2 shows that dissolution changes the porous medium through coupled hydraulic and electrical pathways rather than through porosity alone. Across the three Pe cases, porosity increases from the same initial value, while permeability increases by about 20 times for Pe = 0.1 and by more than 200 times for Pe = 1 and Pe = 10. This permeability growth indicates stronger hydraulic connectivity, but the electrical pathway evolves in the opposite direction for the seismoelectric source strength: pore-fluid concentration and conductivity increase most strongly in the higher Pe cases. Thus, the reactive transport results already contain the central competition of the paper, in which improved flow connectivity can increase hydraulic communication, while enhanced fluid conductivity can screen or reduce electrokinetic coupling.

## Figure 3

Figure 3 translates the reactive transport outputs into the dynamic electrokinetic bridge required by the seismoelectric forward model. As dissolution proceeds, $\sigma(\omega)$ increases with the rise in pore-fluid conductivity, whereas $L(\omega)$ decreases, especially in the Pe = 1 and Pe = 10 cases where concentration increases by roughly 55 to 75 times. This response indicates that the same dissolution process that opens hydraulic pathways also makes the pore fluid more conductive, reducing the relative strength of electrokinetic coupling. The frequency-dependent bridge therefore explains why permeability alone cannot predict the interface EM response: the source term depends on the balance between hydraulic mobility and electrical leakage through $L(\omega)$ and $\sigma(\omega)$.

## Figure 4

Figure 4 shows how the dynamic property changes are transferred into the Schakel-type interface coefficients. The magnitude of $R_E$ decreases to about 0.16 of its initial value for Pe = 0.1, but falls by nearly four orders of magnitude for Pe = 1 and Pe = 10. $T_{\mathrm{TM}}$ shows a similar suppression in the higher Pe cases. These trends indicate that the interface conversion is controlled less by the absolute increase in permeability than by the conductivity-driven reduction of electrokinetic contrast. The phase changes further suggest that dissolution modifies not only the amplitude of the interface response but also the balance between the reflected and transmitted electromagnetic components.

## Figure 5

Figure 5 carries the coefficient-level changes into finite-offset waveform time. The waveform panels are evaluated at selected dissolution states, so dissolution time indexes the evolving material state, while waveform time records the microsecond-scale response after acoustic excitation. The post-$T_0$ interface response weakens strongly as dissolution proceeds, with late-time peak amplitudes decreasing from about $1.6\times10^{15}$ in Pe = 0.1 to about $1.8\times10^{12}$ and $1.3\times10^{12}$ in Pe = 1 and Pe = 10. A frequency-sampling diagnostic controls the pre-$T_0$ residual to about 1.1-1.25% of the post-$T_0$ peak, indicating that the main plotted response is not dominated by pre-arrival leakage. Under these simulations, the waveform result supports the same mechanism inferred from the coefficients: conductivity growth can outweigh the effect of permeability enhancement and substantially reduce the interface EM amplitude.

## Figure 6

Figure 6 examines whether the finite-offset waveforms have the spatial organization expected for an interface EM response. The modeled peak amplitudes vary systematically with receiver position and show polarity changes between the reflected and transmitted sides of the interface. The Liu-type dipole curve is used only as a normalized interpretation of spatial directivity, not as an additional multiplier on the forward waveform. This comparison suggests that the spatial distribution is controlled by both the finite-offset geometry and the evolving interface conversion strength. The polarity behavior is therefore treated as a diagnostic of the modeled source geometry, while the absolute amplitude remains controlled by the dissolution-dependent electrokinetic parameters.

## Figure 7

Figure 7 separates the main parameter pathways using one-at-a-time perturbations. In all three Pe cases, increasing tortuosity-related effects alone produces a small positive contribution, whereas porosity, permeability, and fluid chemistry pathways reduce the waveform peak in the tested configuration. The largest negative contribution comes from fluid chemistry in the Pe = 1 and Pe = 10 cases, with log10 peak changes of about -2.78 and -2.91, compared with full observed changes of about -3.86 and -3.99. This decomposition indicates that the late-stage amplitude loss is not simply a permeability effect. Instead, the strong increase in pore-fluid conductivity and the associated reduction in electrokinetic coupling dominate the simulated response, with nonlinear interactions adding a smaller compensating contribution.

## Figure 8

Figure 8 evaluates normalized response metrics as cautious monitoring indicators within the forward simulations. The normalized $R_E$ metric decreases to about 0.16 for Pe = 0.1 and to about $10^{-4}$ for Pe = 1 and Pe = 10, while the $|L|/|\sigma|$ index falls even more strongly in the higher Pe cases. The high-over-low frequency ratio changes more modestly, generally remaining within about 0.97-1.39 of its initial value for Pe = 1 and about 1.00-1.23 for Pe = 10. These results suggest that normalized amplitude and coupling indices are more sensitive to the conductivity-controlled suppression than the simple spectral ratio used here. They should be interpreted as simulation-based diagnostic quantities rather than evidence of field detectability or experimental validation.

## Chinese Note

这组段落按连续机制链组织，而不是逐图复述。措辞使用 "show", "indicate", "suggest", "under these simulations" 等较克制表达，避免宣称真实实验验证或现场可探测性。Figure 5 已写入频率采样诊断结果：$T_0$ 前 residual 约为 post-$T_0$ 峰值的 1.1-1.25%。
