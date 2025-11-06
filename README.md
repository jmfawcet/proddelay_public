# Pupillometric Production Effect Data

This repository contains the behavioral and pupillometric data associated with:
**Fawcett, J. M., Roberts, B. R. T., Hu, S., Thoms, C., Whitridge, J., MacLeod, C. M., & Willoughby, H.**
*Preparing to produce (without production) is sufficient to elicit a behavioral and pupillometric production effect.* 
_Preprint available at [https://doi.org/10.31234/osf.io/ay6u7_v2](https://doi.org/10.31234/osf.io/ay6u7_v2)_
(Currently under revision at *Psychonomic Bulletin & Review*)


---

## Study Overview

The production effect—better memory for words read aloud than silently—is often attributed to distinctive encoding (e.g., sensorimotor features of speaking). Here, we isolate the role of preparatory attention using a delayed-production paradigm with “Go” signals and catch trials where no response is made.

Across two experiments, participants studied words to be read aloud or silently while withholding responses until a "Go"" signal. On catch trials, the "Go" signal never appeared, letting us test whether preparation to speak (without speaking) enhances memory. Pupillometry (Experiment 2) indexed attention/processing load throughout the trial.

---

## Key Findings

-Behavioral production effect (aloud > silent) replicated across delays; the effect was smaller but present on catch trials, indicating preparation alone confers a memory benefit.

-Pupillometric production effect: pupil dilation was larger for aloud than silent trials, with differences emerging before overt production; a smaller prep-related deflection was also present on catch trials.

-Process evidence: for "Go"" trials, the production effect involved both recollection and familiarity; for catch trials, it was driven primarily by recollection, consistent with a preparatory attention component.

---

## Data Description

The dataset includes:

- **Behavioral data**  
  - Item-level recognition judgments (hits, false alarms).  
  - Confidence ratings
  - Derived sensitivity measures (d′, familiarity, recollection estimates).  

- **Pupillometric data**  
  - Raw pupil size traces sampled via SR Research EyeLink 1000 Plus.  
  - Preprocessed signals (blink removal, interpolation, filtering, baseline correction, z-score normalization).  
  - Trial-wise averaged pupil dilation values (by condition).  

Full data, preprocessing scripts, and analysis code (a zip of the present directory) are also hosted on the **Open Science Framework (OSF): [https://osf.io/dj3nu/](https://osf.io/dj3nu/)**. Note that owing to file size, .edf files are provided only on OSF. These files will need to be added to the relevant data directories or the code used to read the files modified to check for .rds files to gather file names.

---

## Citation

If you use these data, please cite:

> Fawcett, J. M., Roberts, B. R. T., Hu, S., Thoms, C., Whitridge, J., MacLeod, C. M., & Willoughby, H. (2025, September 15). 
> Preparing to produce (without production) is sufficient to elicit a behavioral and pupillometric production effect. 
> *Psychonomic Bulletin & Review* (in revision) Preprint: [https://doi.org/10.31234/osf.io/ay6u7_v2](https://doi.org/10.31234/osf.io/ay6u7_v2)

A machine-readable citation file (`CITATION.cff`) is also included for GitHub integration.

---

## License

[![CC BY-NC 4.0][cc-by-nc-shield]][cc-by-nc]

This work is licensed under a [Creative Commons Attribution-NonCommercial 4.0 International License][cc-by-nc].

[cc-by-nc]: https://creativecommons.org/licenses/by-nc/4.0/
[cc-by-nc-shield]: https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg

The dataset is licensed under **CC BY-NC 4.0**. You are free to use, share, and adapt the materials provided that you give appropriate credit by citing the paper above.

---

## Contact

For questions about the dataset or analysis pipeline, please contact:  
**Jonathan M. Fawcett** – [Memorial University of Newfoundland]

If you use any of this analysis code in your own models, in whole or in part, please cite the above paper. If you find errors, please let me know.

---
