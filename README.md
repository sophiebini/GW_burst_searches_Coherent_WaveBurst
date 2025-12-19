# GW_burst_searches_Coherent_WaveBurst

Coherent WaveBurst ([cWB](https://gwburst.gitlab.io/documentation/latest/html/index.html)) is an open sofware used to search for gravitational-wave (GW) transients without assuming a waveform model [1,2]

This repository contains codes and documentations I've written to search for GW bursts with this algorithm.

Here is a summary what the folders included:

- **autoencoder_for_glitch** to improve the cWB detection efficiency I've developed an autoencoder neural network that identify a know class of noise in GW data (blip glitches) and reject them [3]. These scripts part of the [cWB gitlab project](https://gitlab.com/gwburst) 
- **XGBoost_tutorial**  to improve the cWB detection efficiency a machine learning algorithm, we use a decision tree classifier (XGBoost) to distingish between potential GW signals and noise [4,5]. This notebook shows how to train the classifier using cWB triggers, and apply it to cWB.






## Main references 
1. Sergey Klimenko et al. “A coherent method for detection of gravitational wave bursts”. In: Classical and Quantum Gravity 25.11 (2008), p. 114029.
2. Marco Drago et al. “Coherent WaveBurst, a pipeline for unmodeled gravitational-wave data analysis”. In: SoftwareX 14 (2021), p. 100678.
3. Sophie Bini et al. “An autoencoder neural network integrated into gravitational-wave burst searches to improve the rejection of noise transients”. In: Classical and Quantum Gravity (2023).
4. Tanmaya Mishra et al. “Optimization of model independent gravitational wave search for binary black hole mergers using machine learning”. In: Physical Review D 104.2 (2021),
p. 023014.
5. Marek J. Szczepańczyk et al (with Sophie Bini). “Search for gravitational-wave bursts in the third Advanced LIGO-Virgo run with coherent WaveBurst enhanced by machine learning”. In: Phys.
Rev. D 107 (6 2023), p. 062002. doi: 10 . 1103 / PhysRevD . 107 . 062002
