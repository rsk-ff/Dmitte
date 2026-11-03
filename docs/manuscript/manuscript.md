# Release timing controls the ingestion dose from acute HTO releases to potato and wheat: a crop-growth-driven compartment model with benchmark validation and global sensitivity analysis

*First draft, 2026-11-03. Target journal: Journal of Environmental Radioactivity. Items marked* **[VERIFY]** *must be checked before submission; see `README.md` in this folder.*

Authors: **[to be decided by Fan and supervisor]**

---

## Highlights

- An hourly air–soil–crop tritium model is driven by the WOFOST crop model.
- It is validated against 65 literature benchmarks, with a test set frozen and run once.
- The same 1 h release changes the harvest ingestion dose by about 100 times, depending on its date.
- Constant transfer rates underestimate the worst case by 5–40 times, and the difference comes from crop growth.
- Parameter uncertainty spreads the dose by only 1.2–2 times (Sobol analysis).

(Elsevier limit: 3–5 highlights, at most 85 characters each. **[VERIFY length]**)

---

## Abstract

After an accidental atmospheric release of tritiated water (HTO), the ingestion dose is dominated by organically bound tritium (OBT) in harvested crops. OBT formation depends on photosynthesis and on how dry matter is allocated to the edible organ, both of which change over the growing season. Yet operational assessments often use transfer rates that are constant or that depend only on the conditions at the time of release. We coupled an hourly 11-compartment HTO/OBT model to the WOFOST crop model, so that leaf–air exchange, root water uptake, OBT formation and its allocation follow the simulated crop state, with plant and fruit water balances closed by construction. The model was evaluated against 65 benchmarks compiled from public experiments on potato, wheat and proxy crops. The benchmarks were split by experimental campaign into a calibration set and a test set; the test set was frozen and run once. Three parameters were fitted. Of 25 calibration checks, 24 were within a factor of 2 of the published range, and so were 11 of 13 test checks, including the order-of-magnitude change in grain tritium-loss index with grain-filling stage. A leave-one-group-out cross-validation kept the fitted values stable. For a 1 h release, the harvest dose per kilogram eaten varied by about a factor of 100 (5th–95th percentile) with release date. The same model with season-mean rates varied by only about a factor of 2, underestimating the worst case by a factor of 5 (wheat) to 40 (potato) and overestimating releases before grain or tuber formation by about a factor of 10. Nearly all of this difference came from crop growth rather than the diurnal cycle. In a Sobol analysis of ten parameters, parameter uncertainty spread the dose by only 1.2–2 times, with the OBT isotope discrimination factor dominant for wheat and minimum stomatal resistance for potato. The time of a release within the season therefore matters far more than parameter uncertainty, and dose assessments for acute HTO releases should resolve crop development.

**Keywords:** tritium; organically bound tritium; HTO; crop growth model; WOFOST; accidental release; ingestion dose; Sobol sensitivity analysis

---

## 1. Introduction

Tritium is released to the atmosphere from heavy-water reactors, reprocessing plants and, prospectively, fusion facilities such as ITER and DEMO [refs **[VERIFY]**]. Most atmospheric tritium released in an accident ends up as tritiated water (HTO). It is taken up by crops through leaf exchange within hours, and part of it is incorporated into organic matter as organically bound tritium (OBT). OBT has a longer biological half-life and a higher ingestion dose coefficient than HTO (ICRP, 1995). Because a short release leaves little HTO in plants at harvest, OBT in the edible organ usually dominates the ingestion dose (IAEA, 2014; Galeriu et al., 2013).

Experiments with short HTO exposures have shown that the OBT reaching the harvested organ depends strongly on when the exposure occurs. In wheat, the grain tritium-loss index (TLI) after a 2 h exposure was several times higher during the linear phase of grain filling than before or after it (Diabaté and Strack, 1997; Strack et al., 2005). Similar growth-stage effects have been reported for rice, Chinese cabbage and radish (Choi et al., 2002, 2005) and for vegetables in recent chamber and field studies (Polivkina et al., 2023, 2024 **[VERIFY]**). Time-of-day effects have been shown for soybean and maize exposed by day and by night (Shen and Liu, 2016; Shen and Guo, 2019). The IAEA EMRAS II and MODARIA programmes compiled these experiments and compared models against them (IAEA, 2014, 2022).

Several models link tritium transfer to crop development. UFOTRI and its PLANT-OBT module make OBT formation depend on leaf HTO, plant type and development stage, with night-time formation calibrated separately (Raskob, 1993 **[VERIFY]**; Strack et al., 1995). SOLVEG-II calculates OBT from photosynthesis in a detailed soil–vegetation–atmosphere model (Ota and Nagai, 2011). Melintescu, Galeriu and co-workers used the WOFOST crop model to derive leaf area and photosynthesis for the RODOS tritium module (Melintescu et al., 2002). They then used WOFOST photosynthesis to drive canopy resistance and OBT formation (Melintescu and Galeriu, 2005 **[VERIFY author list]**), and developed this into the CROPTRIT model for winter wheat (Melintescu et al., 2015; Melintescu, 2023). Keum et al. (2005) replaced the constant yield of a compartment model with a time-dependent growth equation and reproduced exposures at different growth stages. Korolevych and Kim (2011, 2014) transferred a tritium module into the CTEM-CLASS land-surface model, which includes plant phenology. Coupling a crop model to OBT formation is therefore not new.

What remains open is how much this coupling matters for dose assessment, and how that compares with the other uncertainties. Reviews of the EMRAS II wheat scenario concluded that the current dynamic models could not yet support a quantitative assessment of uncertainty, and that more tests across crops and conditions were needed (Melintescu and Galeriu, 2017; Galeriu and Melintescu, 2015). Galeriu et al. (1995) sampled UFOTRI parameters and found an uncertainty below one order of magnitude, but did not compare it with the effect of release timing. Probabilistic UFOTRI assessments sample a year of weather (Raskob et al., 2005; Murphy, 2024 **[VERIFY]**), which includes release timing implicitly, but do not separate it. Most model evaluations have been single-site comparisons or intercomparisons on a few scenarios (Barry et al., 1999; IAEA, 2014).

This paper addresses three questions:

1. Can a compartment model whose rates are driven hourly by a general crop model reproduce the published experiments on potato and wheat, when tested on data not used to build it?
2. Within the same model structure, how much of the variation in harvest dose comes from crop growth and how much from the diurnal cycle? How large is the error of constant or release-time-frozen transfer rates?
3. How does parameter uncertainty compare with the effect of release timing?

We answer them with an open model (Dmitte) coupled to WOFOST, a benchmark set of 65 literature values with a test set frozen before evaluation, three non-dynamic baselines built from the same model, and a Morris–Sobol sensitivity analysis.

---

## 2. Methods

### 2.1 Model structure

The model has 11 compartments (Fig. 1): atmospheric HT and soil HT (not used in the HTO analyses here), air moisture HTO, three soil water layers, leaf water HTO, leaf assimilate OBT, stem and leaf structural OBT, fruit (grain or tuber) water HTO and fruit OBT. Transfers are first-order with 33 rates that are updated every hour. Within an hour the system is linear with constant coefficients, so the state is advanced exactly with the matrix exponential of the rate matrix. For stiff hours this avoids the errors of explicit solvers: the repository's earlier Runge–Kutta solver erred by up to 100 % in some hours (Supplementary S1 **[to prepare]**).

**Leaf–air exchange.** HTO exchange between leaf water and air moisture is a symmetric diffusive flux through the total resistance (aerodynamic, boundary layer and canopy stomatal), with the liquid–vapour fractionation factor α = 1.1. The leaf-to-air rate is k_bh,a2 = g ρ_sat /(α W_leaf), where g is the total conductance, ρ_sat the saturation vapour density at leaf temperature and W_leaf the leaf water inventory (IAEA, 2022 **[VERIFY equation number]**). The uptake rate uses the same conductance with the air moisture inventory. Daytime stomatal resistance is the crop's minimum resistance multiplied by the light, water, humidity and temperature factors of UFOTRI. At night, and when relative humidity exceeds 98 %, it is RC_NIGHT times the minimum resistance (Section 2.4). The wind profile is extrapolated from the 10 m measurement height.

**Plant water.** The leaf water inventory is the water of green leaves, from WOFOST's leaf dry weight scaled by LAI/LAI_max and a leaf dry-matter fraction of 0.2. Root water uptake equals net transpiration plus the water needed for growth, so the plant and fruit water balances close by construction. Uptake is taken from the three soil layers in fixed shares (0.2/0.4/0.4). Fruit water is renewed by the net flow needed for fruit growth and by exchange with leaf water (half-life 80 h; not identifiable from the benchmarks and included in the sensitivity analysis). Rain enters the soil only.

**OBT formation.** Following the growth form of IAEA (2022) Eq. 43 **[VERIFY]**, HTO from leaf water enters a fresh-assimilate pool at a rate proportional to WOFOST gross assimilation: ID × 0.6 × gross assimilation × leaf-water activity. Here ID = 0.8 is the isotope discrimination factor and 0.6 the combustion-water equivalent of CH₂O. During daytime this flux is distributed in proportion to radiation. At night it is OBT_DARK = 0.3 times the daytime mean. The assimilate pool empties with residence time TAU_ASSIM, divided according to that day's WOFOST dry-matter partitioning: shoot growth goes to structural OBT, storage-organ growth to fruit OBT, and the remainder (respiration and roots) back to leaf water as HTO. Structural OBT does not turn over (TECDOC-1991 Eq. 43 has growth dilution only). OBT concentrations are expressed per litre of combustion water, using 0.556 kg water per kg dry matter.

**Crop model.** WOFOST 7.2 is run through PCSE 5.5.6 (de Wit et al., 2019 **[VERIFY]**) with the public crop parameters Potato_701 and Winter_wheat_102, the EC2 soil and daily weather from the Wageningen Haarweg station (2004–2008). Hourly forcing is generated from the daily data: temperature follows a cosine between the daily minimum and maximum, radiation follows solar elevation, and relative humidity is recomputed from vapour pressure and hourly temperature. The analyses use potato seasons 2004–2008 (emergence 1 May) and winter wheat seasons 2005–2008 (sown 20 October of the previous year).

**Dose.** Adult ingestion dose coefficients are 1.8 × 10⁻¹¹ Sv Bq⁻¹ for HTO and 4.2 × 10⁻¹¹ Sv Bq⁻¹ for OBT (ICRP, 1995). Doses are reported per kilogram of fresh edible product. Annual doses in Section 3.5 assume 100 kg a⁻¹ of the product, all from the affected location, for comparison only (**[VERIFY / replace with IAEA SRS-19 default intake]**).

### 2.2 Benchmarks and validation design

We compiled 65 scalar benchmarks and 20 experiment-level datasets from the open literature (Table 1; Supplementary S2). The sources are mainly IAEA-TECDOC-1738 (2014), IAEA-TECDOC-1991 (2022), the EMRAS Pickering scenario, Diabaté and Strack (1997), Strack et al. (2005), Korolevych and Kim (2013), Balashov et al. (2011) and Shen and Guo (2019). Each benchmark records:
- the experimental setting (chamber, field box, greenhouse or field);
- the tracer;
- the exposure duration;
- the photosynthetic pathway;
- whether the value is a measurement, a fitted rate, a compilation or a model parameter.

Fifteen scenario cards list the exposure conditions available for driving the model with measured conditions.

The benchmarks were split by experimental campaign, so that no campaign contributes to both sets. The test set (15 benchmarks) consists of three campaigns that had not been used in model development: the Karlsruhe wheat chamber experiments, the Karlsruhe potato and bean chamber experiments, and the Pickering chronic release. One test value (B12) had appeared in earlier internal reports and is a known leak. The test set was fixed by a SHA-256 fingerprint checked by the test suite. It was evaluated once, on 29 October 2026, after all model changes, and the model has not been changed since.

Virtual experiments follow the published protocols. Plants are exposed for 1 h (field box) or 2 h (chamber) to constant air-moisture HTO at 11:00 and at 23:00, every week of every season, and then grow in clean air until harvest. In the chronic case the air HTO is constant over the season. The eight wheat field-box benchmarks are scored on runs driven by the measured temperature, humidity and light of each of the 14 exposures. The in-box aerodynamic resistance was not reported, so 15 s m⁻¹ is used for scoring and 8–30 s m⁻¹ is reported as a range. A comparison is "consistent" if the model median over seasons and exposure dates lies within the published range. It is "close" if it lies within a factor of 2 of the range, the agreement criterion of IAEA (2014, §7.9), and "off" otherwise.

### 2.3 Calibration and cross-validation

Three parameters were set on calibration data only:
- **RC_NIGHT** (night stomatal resistance factor; original value 15, no source) was fitted on seven night-time targets from the field box, greenhouse and review benchmarks, giving 4.
- **TAU_ASSIM** was fitted on the shape of the measured leaf OBT time course in the 14 field-box exposures (2 h, 4 h and 1 d relative to the end of exposure), giving 6 h.
- **OBT_DARK** was set to the central value of benchmark B15, which is therefore not an independent check.

All other parameters come from the literature or from stoichiometry (Table 2).

To check that the fitted values do not depend on one group of experiments, we ran a leave-one-group-out cross-validation on the calibration set, with the scheme committed before the run. RC_NIGHT was refitted with each group of night targets left out: field box, greenhouse, review, and each crop as a whole. TAU_ASSIM was refitted with each field-box exposure left out and with each time-of-day class left out.

### 2.4 Release-timing analysis and frozen-rate baselines

Because the model is linear in the air concentration, the harvest response to a 1 h release at any hour of the season follows from one backward (adjoint) sweep of the hourly propagators. The response to any air concentration history is the superposition of these responses; we verified this against forward runs. For each season we computed the harvest dose per kilogram eaten and the fruit OBT for releases at every hour between emergence and maturity.

To separate the effects of crop growth and of the diurnal cycle, we compared the dynamic model with three baselines that use the same model and the same harvest inventories but freeze part of the time variation:
- **growth frozen**: each hour of day takes its season-mean rates, keeping the diurnal cycle;
- **constant**: the season-mean rates throughout;
- **frozen at release**: the rates of the release hour are kept until harvest, as a quasi-static model evaluated at the time of release would do.

### 2.5 Hypothetical accident scenario

Following the hypothetical scenario of the EMRAS tritium working group **[VERIFY: TC14 chapter 9 reference]**, 3.7 × 10¹⁵ Bq of HTO is released in 1 h at an effective height of 20 m. Three cases are considered: Case 1 with stability class A at 2 m s⁻¹, released at 12:00; Case 2 with class D at 5 m s⁻¹, at 12:00; and Case 3 with class F at 2 m s⁻¹, at 00:00. A Gaussian puff model gives the hourly time-integrated air concentration at 1 m height on the plume axis at 1, 3, 10 and 30 km. Divided by the absolute humidity, this is the air-moisture boundary condition for the compartment model during the whole passage. The stability class is used only for dispersion; the compartment model uses the hourly station weather and neutral stability throughout. Releases were placed on every day between emergence and maturity of every season.

### 2.6 Sensitivity analysis

Ten parameters were varied over uniform ranges (Table 2). Four outputs were considered for each crop, each as the log₁₀ of the geometric mean over seasons:
- the season-mean dose per kilogram eaten;
- the dose for a release at 12:00 in mid filling;
- the dose for the same release at 00:00;
- the season-mean fruit OBT.

Morris screening (r = 20, 220 runs per crop) removed parameters whose μ* was below 10 % of the largest. Sobol first-order and total indices of the remaining eight parameters were computed with N = 256 (2560 runs per crop) using SALib (Herman and Usher, 2017 **[VERIFY]**).

---

## 3. Results

### 3.1 Behaviour of the growth-driven rates

Fig. 2 shows the daily mean of four rates over one wheat and one potato season:
- **Leaf–air exchange** varies within a factor of about 3, following weather and leaf water.
- **OBT formation** follows gross assimilation and declines with senescence.
- **Allocation of assimilate to the storage organ** is zero before grain or tuber formation, rises by two orders of magnitude at the start of storage-organ growth and falls at maturity.

The season-mean (constant) rates therefore overstate allocation before storage-organ formation and understate it during the linear filling phase. This is the main mechanism behind the release-date dependence in Section 3.4.

In a uniform-activity test (all compartments at the same specific activity, no decay) the plant and fruit water compartments stayed constant to numerical precision, confirming the closed water balance. The soil compartments did not, because soil–air exchange is still asymmetric (Section 4.3).

### 3.2 Calibration set

Of 25 calibration checks, 24 were within a factor of 2 of the published range (Fig. 3; Table 3): potato 7 consistent, 4 close, 0 off; wheat 6, 7, 1.

- **Potato:** chronic tuber OBT/air ratio 0.91 (published 0.72–1.14) and leaf exchange rates 0.58 h⁻¹ by day and 0.23 h⁻¹ by night (0.5–1.1 and 0.19–0.42).
- **Wheat field box:** driven by measured conditions, leaf TFWT at the end of a 1 h daytime exposure was 60 % of air-moisture HTO (54–92 %), the night/day ratio 0.35 (0.17–0.38) and grain TLI 0.30 % (0.18–0.32 %).

The one value off was leaf OBT at the end of a night-time exposure (B23, 0.13 % against 0.33–0.43 %). The shape of the leaf OBT time course was also not reproduced. The model keeps forming OBT for 1–3 h after exposure, reaching 1.6 times the measurement at 4 h, whereas the measured OBT had largely formed by the end of exposure (Section 4.3).

### 3.3 Test set and post-test analyses

The frozen test set had model counterparts for 13 of its 15 benchmarks. Eleven of the 13 were within a factor of 2: 3 consistent, 8 close, 2 off (Fig. 3; Table 3).

**Grain-filling stage.** The test set directly tests the growth-driven transfer: grain TLI after a 2 h chamber exposure was 0.031 % before the linear grain-filling phase (published 0.06–0.2 %), 0.41 % during it (0.4–0.9 %) and 0.069 % after it (0.04–0.2 %) (Fig. 4). The model reproduces the order-of-magnitude rise and fall with filling stage, which comes entirely from WOFOST's storage-organ growth. Potato tuber TLI was 0.25 % (0.2–0.3 %).

**The two values off:**
- *B46 (wheat chamber, high light, leaf TFWT):* 2.1 times low. The frozen check used the generic virtual experiment with site wind speed, whereas the chamber was fan-ventilated.
- *B34 (Pickering chronic root-crop OBT/air):* 2.3 times high. The model agrees with the Chalk River garden data in the calibration set (B18), and the two observed ratios themselves differ by a factor of about 3 because they are normalised over different air-averaging windows.

In a post-test analysis, clearly separate from the frozen test, we drove the chamber experiments with their measured hourly conditions (Diabaté and Strack, 1997, Tables 2–3). Leaf TFWT at the end of the daytime exposure rose from 0.52 to 0.89 times the measurement, and the inferred high-light value was about 0.85 times. Most of the B46 bias was therefore due to the experimental conditions not being driven. Night-time leaf uptake remained 1.5 times high, consistent with the field box.

**Cross-validation.** Refitting RC_NIGHT without each group of targets gave values of 4 or 5. All 14 left-out targets stayed within a factor of 2, with an RMS log₁₀ miss of 0.071 out of sample against 0.047 in sample. TAU_ASSIM stayed at 6 h in all 14 leave-one-exposure-out folds; leaving out the daytime exposures as a class gave 4 h. The large error in the shape of the leaf OBT time course (RMS log₁₀ 0.55) was the same in and out of sample, so it reflects model structure rather than over-fitting.

### 3.4 Release-date dependence and frozen-rate baselines

With the dynamic model, the harvest dose per kilogram eaten from a 1 h release varied by a factor of 100 (potato) and 92 (wheat) between the 5th and 95th percentiles of release hours (Fig. 5; Table 4). With constant rates the factor was only 2.1–2.2. Relative to the dynamic model, the constant baseline:
- underestimated the season-mean dose by a factor of 7 (potato) and 2 (wheat);
- underestimated the season maximum by a factor of 42 (potato) and 5 (wheat);
- overestimated releases before storage-organ formation by a factor of 5–10.

Binned by release date, the dynamic/constant ratio for wheat was 0.1–0.2 before anthesis, 2.4–5 during the first 30 days after anthesis and 1.5–1.7 later. For potato it rose from 0.1–0.5 before tuber initiation to 8–13 late in tuber bulking.

**Growth versus the diurnal cycle.** Almost all of the difference is due to crop growth. The growth-frozen baseline, which keeps the diurnal cycle, differed from the constant baseline by only a factor of 1.0–1.3.

**Rates frozen at release.** Freezing the rates at the time of release was not a remedy:
- it gives zero dose for releases before the grain or tuber exists;
- it overestimates the wheat maximum by a factor of 1.45;
- it produces a noon/midnight dose ratio of 3.4, whereas in the dynamic model the ratio is about 1, because HTO retained in leaves after a night release is fixed by photosynthesis the next morning.

**Late potato releases.** Part of the high dose from late potato releases is HTO in tuber water at harvest. For OBT alone, the season maximum over the median was 5 for potato and 13 for wheat.

### 3.5 Hypothetical accident

At 1 km, the median annual ingestion dose over release dates and seasons was 0.18–4.4 mSv for potato and 0.22–5.4 mSv for wheat, increasing from Case 1 (A) to Case 3 (F) (Fig. 6; Table 5). The factor of about 25 between Case 3 and Case 1 comes from the time-integrated air concentration. Within one case, release date changed the dose by about 70 times (5th–95th percentile), more than the difference between cases. The constant-rate model compressed this range to 1.3–1.6 times and underestimated the worst case by 7–70 times. Distance changes the time-integrated concentration but leaves the ratios between models almost unchanged.

### 3.6 Sensitivity analysis

Morris screening removed the liquid–vapour fractionation factor and structural OBT turnover. Sobol indices of the remaining eight parameters (Fig. 7; Table 6) showed weak interactions (S1 ≈ ST).

- **Wheat:** the OBT isotope discrimination factor ID dominated every output (ST 0.76–0.92).
- **Potato dose:** the minimum stomatal resistance dominated (ST 0.42–0.69). For releases during tuber filling, the leaf-to-fruit water exchange KBH_FH also mattered (ST 0.27–0.30). Potato OBT was controlled by ID (ST 0.97).
- **Calibrated parameters:** RC_NIGHT and TAU_ASSIM had ST ≤ 0.05 for every output, and OBT_DARK at most 0.07.

The 95th/5th percentile ratio of the outputs over the parameter ranges was 1.24–1.95. This is two orders of magnitude smaller than the release-date effect and smaller than the dynamic/constant differences.

---

## 4. Discussion

### 4.1 Release timing as the leading factor

Over the realistic parameter ranges explored here, the date of an acute HTO release within the growing season matters more for the harvest ingestion dose than any model parameter. The mechanism is simple: OBT reaches the edible organ only through allocation of fresh assimilate to that organ, and that allocation is concentrated in a few weeks of linear grain or tuber growth. A model with constant rates spreads this window over the season, underestimating releases in the sensitive window and overestimating releases before it.

This agrees with the experimental finding that grain OBT is highest for exposures during linear grain filling (Diabaté and Strack, 1997; Strack et al., 2005), and with the growth-stage dependence reported by Keum et al. (2005) and Choi et al. (2005). Our contribution is to quantify the consequence for the dose and to attribute it, within one model structure, to crop growth rather than to diurnal variation. For emergency planning, the factor of 5–40 between the worst case and a constant-rate estimate is larger than the factor of 2 usually accepted as model–data agreement for this class of models (IAEA, 2014).

These results also bear on the choice between operational and research models. UFOTRI and CROPTRIT already include development-stage effects, but through crop-specific parameters or a site-adapted growth model (Raskob, 1993 **[VERIFY]**; Melintescu, 2023). Our results suggest that the structural element that matters most is the timing of storage-organ growth. A model that represents this timing, even with modest physiological detail, captures most of the release-date variation. Conversely, quasi-static approaches that freeze conditions at the time of release fail both before storage-organ formation and in the day–night contrast.

### 4.2 Parameter uncertainty and calibration

The Sobol analysis agrees with Galeriu et al. (1995), who found parameter uncertainty in UFOTRI below one order of magnitude with stomatal resistance and LAI among the leading parameters. It adds a direct comparison with release timing. That the calibrated parameters (RC_NIGHT, TAU_ASSIM) have negligible influence on the dose means their calibration errors do not propagate to the assessment endpoint. The isotope discrimination factor, taken from the literature, dominates wheat outputs, and better constraining it would reduce dose uncertainty most.

Two caveats apply:
- The sensitivity analysis covers parameter values, not the functional forms. The chamber analysis suggests that the night-time formation of OBT cannot be described by a constant fraction of daytime formation (Section 4.3).
- The ranges are literature or plausibility ranges rather than calibrated posteriors.

### 4.3 Limitations

- **Leaf OBT short-term dynamics.** The single first-order assimilate pool cannot reproduce both the rapid labelling at the end of exposure and the 1 d decline (B22, B23). The measurements suggest labelling faster than gross photosynthesis, for example of exchangeable hydrogen in sugars. Satoh et al. (2025 **[VERIFY]**) and Galeriu and Melintescu (2017) point to non-photosynthetic and night-time processes. This affects leaf OBT shortly after exposure more than OBT in the harvested organ, for which the TLI benchmarks agree.
- **Night-time and dawn/dusk stomatal response.** Night-time leaf uptake is about 1.5 times high in both the field box and the chamber, so night-release doses (Case 3) may be overestimated. The model also responds only to instantaneous light, so stomata open too little at dawn and too much at dusk.
- **Soil water.** Soil–air exchange is still asymmetric (deposition via deposition velocity, re-emission via net evaporation), and the soil compartments do not pass the uniform-activity test. Soil contributes little to harvest OBT after a short release but more for chronic releases.
- **Site, varieties and forcing.** One site (Wageningen) and one variety per crop were used, with hourly forcing synthesised from daily weather. The chamber experiments used spring wheat, while WOFOST simulates winter wheat.
- **Atmospheric coupling.** Air is a boundary condition: re-emitted HTO does not return to air, and the dispersion model uses the case wind speed while leaf exchange uses the station wind speed. The dispersion part was not validated.
- **Benchmark coverage.** Some benchmarks are from proxy crops, and several published ranges are narrow relative to experimental error.

### 4.4 Reproducibility

All code, the benchmark set, scenario cards, the frozen-split fingerprint and the scripts producing every number and figure are public (Section "Data availability"). All inputs are downloaded from fixed public sources by one script. The full validation and analysis run in about one hour on a laptop.

---

## 5. Conclusions

- An hourly compartment model whose transfer rates are driven by WOFOST reproduced 24 of 25 calibration and 11 of 13 frozen test benchmarks within a factor of 2, including the order-of-magnitude dependence of grain TLI on grain-filling stage.
- The harvest ingestion dose from a 1 h HTO release varied by about a factor of 100 with release date. A constant-rate version of the same model underestimated the worst case by 5–40 times and overestimated early releases by about 10 times. The effect is due to crop growth, not to the diurnal cycle.
- Parameter uncertainty changed the dose by only 1.2–2 times, so release timing is the leading factor in dose assessment for acute HTO releases.
- Operational models for accidental tritium releases should represent the timing of storage-organ growth.

---

## Data availability

Code, benchmarks and results: https://github.com/rsk-ff/Dmitte (archived at Zenodo, DOI **[to be created]**). Weather, soil and crop parameter inputs are public and downloaded by `validation/build_inputs.py`.

## Declaration of generative AI use

During the preparation of this work the authors used Claude (Anthropic) to assist with code development, analysis and drafting of the manuscript. After using this tool, the authors reviewed and edited the content as needed and take full responsibility for the content of the publication. **[VERIFY wording against the current Elsevier policy]**

## CRediT authorship contribution statement / Funding / Declaration of competing interest

**[to be completed]**

---

## Tables

**Table 1.** Benchmark set by crop, process and set. **[to generate from `validation/benchmarks/literature_benchmarks.csv` and `split.csv`]**

**Table 2.** Model parameters, values, sources and sensitivity ranges.

| Parameter | Value | Source | SA range |
|---|---|---|---|
| α, liquid–vapour HTO fractionation | 1.1 | literature | 1.08–1.12 |
| Leaf dry-matter fraction | 0.2 | IAEA (2022); Diabaté and Strack (1997) | 0.15–0.25 |
| ID, isotope discrimination of OBT formation | 0.8 | IAEA (2022) Eq. 43 | 0.6–1.0 |
| Combustion water of dry matter | 0.556 kg kg⁻¹ | stoichiometry (6.2 % H) | — |
| Combustion water of CH₂O | 0.6 kg kg⁻¹ | stoichiometry | — |
| RC_NIGHT, night stomatal resistance factor | 4 | calibrated (§2.3) | 3–6 |
| TAU_ASSIM, assimilate residence time | 6 h | calibrated (§2.3) | 2–12 h |
| OBT_DARK, night/day OBT formation | 0.3 | B15 central value; IAEA (2014) 0.1–0.5 | 0.1–0.5 |
| Leaf→fruit water exchange half-life | 80 h | not identifiable | 20–200 h |
| Structural OBT turnover | 0 | IAEA (2022) Eq. 43 | 0–ln2/10 d |
| Minimum stomatal resistance | crop value | UFOTRI | ×0.5–1.5 |
| Soil percolation scheme | constant flux | — | constant / Darcy |
| Dose coefficients HTO / OBT | 1.8 / 4.2 × 10⁻¹¹ Sv Bq⁻¹ | ICRP (1995) | — |

**Table 3.** Calibration and test results by benchmark. **[from `validation/results/benchmark_comparison*.csv`]**

**Table 4.** Season statistics of the harvest dose per kg eaten for the dynamic model and the three baselines (median over seasons; mSv kg⁻¹ per Bq h kg⁻¹).

| Crop | Model | Season mean | Season max | p95/p5 | Noon/midnight | Mean / dynamic | Max / dynamic |
|---|---|---|---|---|---|---|---|
| Potato | dynamic | 9.9e-12 | 8.9e-11 | 100 | 0.97 | 1 | 1 |
| | growth frozen | 1.5e-12 | 2.7e-12 | 2.7 | 1.47 | 0.15 | 0.030 |
| | constant | 1.4e-12 | 2.1e-12 | 2.1 | 0.99 | 0.14 | 0.024 |
| | frozen at release | 6.1e-12 | 4.7e-11 | — | 1.13 | 0.61 | 0.54 |
| Wheat | dynamic | 8.2e-12 | 3.6e-11 | 92 | 1.02 | 1 | 1 |
| | growth frozen | 4.3e-12 | 8.5e-12 | 3.0 | 1.51 | 0.53 | 0.24 |
| | constant | 4.1e-12 | 6.6e-12 | 2.2 | 1.00 | 0.50 | 0.19 |
| | frozen at release | 6.7e-12 | 5.2e-11 | — | 3.4 | 0.82 | 1.45 |

**Table 5.** Annual ingestion dose at 1 km in the hypothetical accident (mSv; 5 % / median / 95 % / max over release dates and seasons).

| Crop | Case | Dynamic | Constant | Frozen at release |
|---|---|---|---|---|
| Potato | 1 (A, 12:00) | 0.014 / 0.18 / 0.97 / 3.9 | 0.043 / 0.049 / 0.056 / 0.058 | 0 / 0.13 / 0.63 / 1.7 |
| | 2 (D, 12:00) | 0.071 / 0.91 / 4.9 / 19 | 0.22 / 0.25 / 0.28 / 0.29 | 0 / 0.66 / 3.2 / 8.6 |
| | 3 (F, 00:00) | 0.33 / 4.4 / 29 / 99 | 1.0 / 1.2 / 1.4 / 1.4 | 0 / 3.0 / 17 / 36 |
| Wheat | 1 (A, 12:00) | 0.013 / 0.22 / 1.0 / 1.7 | 0.13 / 0.15 / 0.20 / 0.20 | 0 / 0.040 / 1.3 / 2.8 |
| | 2 (D, 12:00) | 0.064 / 1.1 / 5.1 / 8.6 | 0.65 / 0.77 / 1.0 / 1.0 | 0 / 0.20 / 6.7 / 14 |
| | 3 (F, 00:00) | 0.38 / 5.4 / 18 / 31 | 3.1 / 3.7 / 4.9 / 4.9 | 0 / 0.38 / 9.8 / 15 |

**Table 6.** Sobol total indices.

| Parameter | Potato dose mean | fill noon | fill night | OBT mean | Wheat dose mean | fill noon | fill night | OBT mean |
|---|---|---|---|---|---|---|---|---|
| RCMIN factor | 0.69 | 0.42 | 0.55 | 0.02 | 0.06 | 0.00 | 0.07 | 0.04 |
| ID | 0.15 | 0.19 | 0.07 | 0.97 | 0.89 | 0.90 | 0.76 | 0.92 |
| KBH_FH | 0.05 | 0.27 | 0.30 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| Percolation | 0.06 | 0.01 | 0.03 | 0.00 | 0.04 | 0.00 | 0.00 | 0.02 |
| OBT_DARK | 0.00 | 0.01 | 0.01 | 0.00 | 0.00 | 0.06 | 0.07 | 0.00 |
| TAU_ASSIM | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.05 | 0.00 |
| Leaf dry-matter fraction | 0.00 | 0.05 | 0.01 | 0.00 | 0.00 | 0.03 | 0.04 | 0.00 |
| RC_NIGHT | 0.01 | 0.01 | 0.00 | 0.00 | 0.00 | 0.00 | 0.02 | 0.00 |
| Output p95/p5 | 1.24 | 1.71 | 1.95 | 1.59 | 1.52 | 1.62 | 1.68 | 1.62 |

## Figure captions

**Fig. 1.** Model compartments and main HTO/OBT pathways. Rates of the blue compartments are updated hourly from WOFOST. (`figures/fig1_structure.png`)

**Fig. 2.** Daily mean of four growth-driven transfer rates over one wheat and one potato season (Wageningen 2006). Dashed lines show the season means used by the constant-rate baseline; the grey band is the linear phase of storage-organ growth. (`figures/fig2_rates.png`)

**Fig. 3.** Model–benchmark agreement for the 25 calibration and 13 test checks: signed log₁₀ of the factor by which the model median misses the published range (0 inside the range). The grey band marks a factor of 2. (`figures/fig3_benchmarks.png`)

**Fig. 4.** Grain TLI after a 2 h exposure as a function of exposure date (left) and by grain-filling phase against the test benchmarks B53–B55 (right; grey boxes are the published ranges). (`figures/fig4_tli_phase.png`)

**Fig. 5.** Harvest dose per kilogram eaten for a 1 h release at each date of the season: dynamic model and frozen-rate baselines. (`../../analysis/results/figures/fig_b1_release_dates.png`)

**Fig. 6.** Annual ingestion dose at 1 km in the hypothetical accident by release date, for the three dispersion cases. (`../../analysis/results/figures/fig_b2_emras.png`)

**Fig. 7.** Morris screening and Sobol total indices. (`../../analysis/results/figures/fig_b3_sensitivity.png`)

---

## References

**[VERIFY all entries: volumes, pages, DOIs, and author lists; complete the ones marked]**

- Balashov, Yu.S., et al., 2011. *Fusion Sci. Technol.* 60, 1215–1219.
- Barry, P.J., Watkins, B.M., Belot, Y., Davis, P.A., Edlund, O., Galeriu, D., Raskob, W., et al., 1999. Intercomparison of model predictions of tritium concentrations in soil and foods following acute airborne HTO exposure. *J. Environ. Radioact.* 42, 191–207. doi:10.1016/S0265-931X(98)00054-X
- Choi, Y.H., et al., 2002. TFWT and OBT in the rice plant acutely exposed to atmospheric HTO vapor. *J. Environ. Radioact.*
- Choi, Y.H., et al., 2005. Tritium levels in Chinese cabbage and radish plants acutely exposed to HTO vapor at different growth stages. *J. Environ. Radioact.* doi:10.1016/j.jenvrad.2005.04.004
- de Wit, A., et al., 2019. 25 years of the WOFOST cropping systems model. *Agric. Syst.* 168, 154–167.
- Diabaté, S., Strack, S., 1997. Organically bound tritium in wheat after short-term exposure to atmospheric tritium under laboratory conditions. *J. Environ. Radioact.* 36, 157–175.
- Galeriu, D., Davis, P., Chouhan, S., Raskob, W., 1995. Uncertainty and sensitivity analysis for the environmental tritium code UFOTRI. *Fusion Technol.* 28. doi:10.13182/FST95-A30511
- Galeriu, D., Melintescu, A., 2015. Progresses in tritium accident modelling in the frame of IAEA EMRAS II. *Fusion Sci. Technol.* 67. doi:10.13182/FST14-T26
- Galeriu, D., Melintescu, A., 2017. Relevance of night production of OBT in crops. *Fusion Sci. Technol.* 71, 595–599. doi:10.1080/15361055.2016.1273705
- Galeriu, D., Melintescu, A., Strack, S., Atarashi-Andoh, M., Kim, S.B., 2013. An overview of organically bound tritium experiments in plants following a short atmospheric HTO exposure. *J. Environ. Radioact.* 118, 40–56. doi:10.1016/j.jenvrad.2012.11.005
- Herman, J., Usher, W., 2017. SALib: an open-source Python library for sensitivity analysis. *J. Open Source Softw.* 2, 97.
- IAEA, 2014. Transfer of Tritium in the Environment after Accidental Releases from Nuclear Facilities. IAEA-TECDOC-1738. IAEA, Vienna.
- IAEA, 2022. Harmonization and Intercomparison of Models for Accidental Tritium Releases to the Atmosphere. IAEA-TECDOC-1991. IAEA, Vienna.
- ICRP, 1995. Age-dependent doses to members of the public from intake of radionuclides: Part 5. ICRP Publication 72. *Ann. ICRP* 26(1).
- Keum, D.K., Lee, H.S., Kang, H.S., Jun, I., Choi, Y.H., Lee, C.W., 2005. Prediction of tritium level in agricultural plants after short term exposure to HTO vapor and its comparison with experimental results. *Health Phys.* doi:10.1097/01.hp.0000174811.24386.cc
- Korolevych, V., Kim, S.B., 2011. Modelling and validation of OBT formation in tomato and potato plants. *Fusion Sci. Technol.* 60.
- Korolevych, V.Y., Kim, S.B., 2013. *J. Environ. Radioact.* 118, 113–120.
- Korolevych, V., Kim, S.B., Davis, P., 2014. OBT/HTO ratio in agricultural produce subject to routine atmospheric releases of tritium. *J. Environ. Radioact.*
- Melintescu, A., 2023. Radiological impact assessment of acute tritium releases — a soil dynamic model. *Fusion Sci. Technol.* doi:10.1080/15361055.2023.2223723
- Melintescu, A., Galeriu, D., 2005 (author list to verify: Barescut, J.C., Gariel, J.C., Péres, J.M., Melintescu, A., Galeriu, D.). A versatile model for tritium transfer from atmosphere to plant and soil. *Radioprotection* 40(S1), S437–S442. doi:10.1051/radiopro:2005s1-064
- Melintescu, A., Galeriu, D., 2017. Uncertainty of current understanding regarding OBT formation in plants. *J. Environ. Radioact.* doi:10.1016/j.jenvrad.2016.11.026
- Melintescu, A., Galeriu, D., Diabaté, S., Strack, S., 2015. Preparatory steps for a robust dynamic model for organically bound tritium dynamics in agricultural crops. *Fusion Sci. Technol.* 67, 479–482. doi:10.13182/FST14-T59
- Melintescu, A., Galeriu, D., Marica, A., 2002. Using WOFOST crop model for data base derivation of tritium and terrestrial food chain modules in RODOS. *Radioprotection* 37. doi:10.1051/radiopro/2002154
- Murphy, J., 2024. Evaluating postulated tritium releases on a large DOE site using UFOTRI. OSTI report.
- Ota, M., Nagai, H., 2011. Development and validation of a dynamical atmosphere–vegetation–soil HTO transport and OBT formation model. *J. Environ. Radioact.* 102, 813–823.
- Polivkina, E., et al., 2024. Tritium uptake in crops in the area with a high level of atmospheric tritium oxide in the territory of the former Semipalatinsk test site. *PLoS ONE*.
- Raskob, W., 1993. Description of the tritium model UFOTRI. KfK report **[VERIFY]**.
- Raskob, W., Velarde, M., Perlado, J.M., 2005. Tritium dose assessments for the potential European ITER site Vandellós. *Fusion Sci. Technol.*
- Satoh, Y., et al., 2025. Necessity of considering non-photosynthetic metabolism to reproduce the behavior of OBT in a leafy vegetable. *J. Hazard. Mater.* doi:10.1016/j.jhazmat.2025.139518
- Shen, H., Guo, X., 2019. *Appl. Radiat. Isot.* 154, 108905.
- Shen, H., Liu, X., 2016. Tritium concentration in soybean plants exposed to atmospheric HTO during nighttime and daytime. *Nucl. Sci. Tech.*
- Strack, S., Diabaté, S., Müller, H., Raskob, W., 1995. Organically bound tritium formation and translocation in crop plants: modelling and experimental results. *Fusion Technol.* 28.
- Strack, S., Diabaté, S., Raskob, W., 2005. *Fusion Sci. Technol.* 48, 767–770.
