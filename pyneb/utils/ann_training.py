"""

Save a simple ANN 


@authors: ChristopheMorisset & RogelioOrozcoDuarte & OskarArangure

Routine flow:
1. Crear un ion de PyNeb
2. Generar puntos aleatorios de Te y ne
3. Calcular el ratio diagnóstico con getEmissivity()
4. Construir X_train y y_train
5. Entrenar la ANN usando manage_RM
6. Guardar la ANN
7. Guardar metadata en JSON <-- this is to document how the ANN was trained.
"""

import os
import json
import numpy as np
import pyneb as pn
from pyneb.utils.ai_neb import manage_RM


def train_temden_ann(atom, ion, to_eval, mode='tem',
                     filename='ANN_model',
                     n_train=1000,
                     tem_min=3000.,
                     tem_max=20000.,
                     den_min=10.,
                     den_max=1e6,
                     scaling=True,
                     use_log=True,
                     random_seed=0,
                     save_train=True,
                     save_metadata=True):
    """
    Train and save an ANN model for PyNeb Te/ne diagnostics.

    Parameters
    ----------
    atom, ion : str, int
        Ion definition, e.g. atom='O', ion=3.

    to_eval : str
        PyNeb expression for the diagnostic ratio,
        e.g. 'L(4363)/L(5007)'.

    mode : {'tem', 'den'}
        If 'tem', train ANN to predict Te from ratio and ne.
        If 'den', train ANN to predict ne from ratio and Te.

    filename : str
        Output filename without extension.

    n_train : int
        Number of training points.

    tem_min, tem_max : float
        Temperature range in K.

    den_min, den_max : float
        Density range in cm^-3.

    scaling : bool
        Whether to use input scaling in manage_RM.

    use_log : bool
        Whether manage_RM applies log scaling internally.

    random_seed : int
        Random seed for reproducibility.

    save_train : bool
        Whether to save the training data in the ANN file.

    save_metadata : bool
        Whether to save a JSON metadata file next to the ANN file.

    Returns
    -------
    RM : manage_RM
        Trained regression model manager.
    """

#raises an error if we dont select an expected mode
    if mode not in ('tem', 'den'):
        raise ValueError("mode must be 'tem' or 'den'")

    A = pn.Atom(atom, ion)
    rng = np.random.default_rng(random_seed) #generates random points with a fixed seed

    # Training grid
    tem_train = tem_min + (tem_max - tem_min) * rng.random(n_train)

    den_train = 10**(
        np.log10(den_min)
        + (np.log10(den_max) - np.log10(den_min)) * rng.random(n_train)
    )

    #get the emissivity for our Atom
    def L(wave):
        return A.getEmissivity(
            tem=tem_train,
            den=den_train,
            wave=wave,
            product=False
        )
    #get the levels
    def I(lev_i, lev_j):
        return A.getEmissivity(
            tem=tem_train,
            den=den_train,
            lev_i=lev_i,
            lev_j=lev_j,
            product=False
        )

    safe_namespace = {
        "L": L,
        "I": I,
        "np": np
    }

    ratio_train = eval(
        to_eval,
        {"__builtins__": {}},
        safe_namespace
    )

    ratio_train = np.asarray(ratio_train)

    valid = np.isfinite(ratio_train) & (ratio_train > 0)
    #remeber y (temperature/density) is the output, x is our ratio + (temperature/density)
    if mode == 'tem':
        X = np.asarray((ratio_train[valid], den_train[valid])).T
        y = np.log10(tem_train[valid])
    else:
        X = np.asarray((ratio_train[valid], tem_train[valid])).T
        y = np.log10(den_train[valid])

    #call to manage_RM
    RM = manage_RM(
        X_train=X,
        y_train=y,
        scaling=scaling,
        use_log=use_log,
        random_seed=random_seed
    )

    RM.init_RM(
        solver='lbfgs',
        activation='tanh',
        hidden_layer_sizes=(10, 10),
        tol=1e-6,
        max_iter=20000
    )
    #training
    RM.train_RM()
    #save the ANN
    RM.save_RM(filename, save_train=save_train)

    if save_metadata:
        metadata = {
            "atom": atom,
            "ion": ion,
            "to_eval": to_eval,
            "mode": mode,
            "tem_min": tem_min,
            "tem_max": tem_max,
            "den_min": den_min,
            "den_max": den_max,
            "n_train_requested": n_train,
            "n_train_valid": int(valid.sum()),
            "scaling": scaling,
            "use_log": use_log,
            "random_seed": random_seed,
            "pyneb_version": getattr(pn, "__version__", None),
            "filename": filename + ".ai4neb_sk"
        }

        with open(filename + ".json", "w") as f:
            json.dump(metadata, f, indent=2)

    return RM



 #routine to train several classic diagnostics
import os


def train_default_temden_anns(output_dir='.',
                              diagnostics=None,
                              n_train=1000,
                              tem_min=5000.,
                              tem_max=20000.,
                              den_min=10.,
                              den_max=1e5,
                              scaling=True,
                              use_log=True,
                              random_seed=0,
                              save_train=True,
                              save_metadata=True,
                              overwrite=True):
    """
    Train and save a set of ANN models for PyNeb temperature and
    density diagnostics.

    If diagnostics is None, a default set of common diagnostics is used.
    Otherwise, diagnostics must be a list of dictionaries.

    Each diagnostic dictionary must contain:

        name    : output filename without extension
        atom    : atomic symbol, e.g. 'O'
        ion     : ionization stage, e.g. 3
        to_eval : PyNeb diagnostic expression, e.g. 'L(4363)/L(5007)'
        mode    : 'tem' or 'den'

    Parameters
    ----------
    output_dir : str
        Directory where the ANN files will be saved.

    diagnostics : list of dict or None
        List of diagnostics to train. If None, use the default diagnostics.

    n_train : int
        Number of training points for each ANN.

    tem_min, tem_max : float
        Temperature range in K.

    den_min, den_max : float
        Density range in cm^-3.

    scaling : bool
        Whether to use input scaling in manage_RM.

    use_log : bool
        Whether manage_RM applies log scaling internally.

    random_seed : int
        Random seed for reproducibility.

    save_train : bool
        Whether to save the training data in the ANN file.

    save_metadata : bool
        Whether to save JSON metadata next to the ANN file.

    overwrite : bool
        If False, skip models whose .ai4neb_sk file already exists.

    Returns
    -------
    trained_models : dict
        Dictionary containing the trained manage_RM objects.
        Keys are the diagnostic names.
    """

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    default_diagnostics = [
        # Temperature diagnostics
        {
            "name": "O3_4363_5007_tem_ann",
            "atom": "O",
            "ion": 3,
            "to_eval": "L(4363)/L(5007)",
            "mode": "tem"
        },
        {
            "name": "N2_5755_6584_tem_ann",
            "atom": "N",
            "ion": 2,
            "to_eval": "L(5755)/L(6584)",
            "mode": "tem"
        },
        {
            "name": "S3_6312_9069_tem_ann",
            "atom": "S",
            "ion": 3,
            "to_eval": "L(6312)/L(9069)",
            "mode": "tem"
        },

        # Density diagnostics
        {
            "name": "S2_6731_6716_den_ann",
            "atom": "S",
            "ion": 2,
            "to_eval": "L(6731)/L(6716)",
            "mode": "den"
        },
        {
            "name": "O2_3726_3729_den_ann",
            "atom": "O",
            "ion": 2,
            "to_eval": "L(3726)/L(3729)",
            "mode": "den"
        },
        {
            "name": "Cl3_5538_5518_den_ann",
            "atom": "Cl",
            "ion": 3,
            "to_eval": "L(5538)/L(5518)",
            "mode": "den"
        },
        {
            "name": "Ar4_4740_4711_den_ann",
            "atom": "Ar",
            "ion": 4,
            "to_eval": "L(4740)/L(4711)",
            "mode": "den"
        },
    ]

    if diagnostics is None:
        diagnostics = default_diagnostics

    required_keys = ["name", "atom", "ion", "to_eval", "mode"]

    trained_models = {}

    for diag in diagnostics:

        for key in required_keys:
            if key not in diag:
                raise ValueError(
                    "Each diagnostic must contain the key '{}'. "
                    "Problematic diagnostic: {}".format(key, diag)
                )

        if diag["mode"] not in ("tem", "den"):
            raise ValueError(
                "Diagnostic '{}' has invalid mode '{}'. "
                "mode must be 'tem' or 'den'.".format(
                    diag["name"], diag["mode"]
                )
            )

        filename = os.path.join(output_dir, diag["name"])
        ann_file = filename + ".ai4neb_sk"

        if os.path.exists(ann_file) and not overwrite:
            print("Skipping existing ANN model: {}".format(ann_file))
            continue

        print("Training ANN model: {}".format(diag["name"]))
        print("  Ion: {}{}".format(diag["atom"], diag["ion"]))
        print("  Ratio: {}".format(diag["to_eval"]))
        print("  Mode: {}".format(diag["mode"]))

        RM = train_temden_ann(
            atom=diag["atom"],
            ion=diag["ion"],
            to_eval=diag["to_eval"],
            mode=diag["mode"],
            filename=filename,
            n_train=n_train,
            tem_min=tem_min,
            tem_max=tem_max,
            den_min=den_min,
            den_max=den_max,
            scaling=scaling,
            use_log=use_log,
            random_seed=random_seed,
            save_train=save_train,
            save_metadata=save_metadata
        )

        trained_models[diag["name"]] = RM

        print("  Saved: {}".format(ann_file))
        print()

    return trained_models