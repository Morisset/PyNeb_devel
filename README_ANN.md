# ANN Integration in PyNeb

This document describes the integration of the Artificial Neural Network (ANN)
regression utilities into PyNeb. The goal is to make the ANN functionality
self-contained inside PyNeb without requiring the external `ai4neb` package.

## Overview

Previously, PyNeb relied on the external package **ai4neb** to provide the
`manage_RM` class used for ANN-based regression models. This dependency has
been removed by integrating a simplified and adapted version of `manage_RM`
directly into PyNeb.

The ANN utilities are now located in: pyneb/utils/ai_neb.py


This allow Pyneb to perform ANN-based diagnostic without requiering external packages.

---

# Files modified

The following files were modified to integrate the aNN functionality:

### 1. 'pyneb/utils/ai_neb.py' 

Contains the ANN regression manage clase:
manage_RM

This class handles:

- training neural network regressors
- prediction
- scaling of input data
- saving and loading trained models

Only the **scikit-learn ANN implementation ("SK_ANN")** is currently supported

--

### 2/ 'pyneb/core/diags.py'

ANN-based diagnostic use the integrated 'manage_RM' class.

The training workflow is:


1.- Generate synthetic training data (temperature, density)
2.- Compute emissivites using PyNeb atomic models
3.- Train an ANN to approximate the inverse mapping
4.- Use the ANN to predcit temperature or density

The following functio uses the ANN mode:

_getPopulations_ANN()

---

### 3. 'pyneb/core/pynebcore.py'

The ANN solver was integrated into the 'Atom' class via:

_getTemDen_ANN()

This method allows computing temperature or density from line ratios using neural network approximation.

Example:

    python
O3.getTemDen(
    0.01,
    den=1000,
    to_eval='L(4363)/L(5007)',
    method='ANN'
)

Usage Example:

import pyneb as pn

O3 = pn.Atom('O', 3)

tem = O3.getTemDen(
    0.01,
    den=1000,
    to_eval='L(4363)/L(5007)',
    method='ANN'
)

print(tem)
