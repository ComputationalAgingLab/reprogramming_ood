import numpy as np

def _compute_freqs(expected, observed, nbins, bucket_type='bins'):
    if bucket_type == "bins":
        breakpoints = np.histogram(expected, nbins)[1]
    elif bucket_type == "quantiles":
        breakpoints = np.arange(0, nbins + 1) / (nbins) * 100
        breakpoints = np.percentile(expected, breakpoints)
    elif bucket_type == "beta":
        breakpoints = np.histogram([0, 1], nbins)[1]
    else:
        raise NotImplementedError()
    # Calculate frequencies
    expected_freqs = np.histogram(expected, breakpoints)[0] / len(expected)
    observed_freqs = np.histogram(observed, breakpoints)[0] / len(observed)
    return expected_freqs, observed_freqs

def _psi(expected: np.ndarray, 
         observed: np.ndarray, 
         nbins: int = 16,
         bucket_type: str = "bins", 
         a_min=0.0001,
         p_max=1.0
         ) -> float:
    """Calculate PSI metric for two arrays.
    
    Parameters
    ----------
        expected : list-like
            Array of expected values
        actual : list-like
            Array of actual values
        bucket_type : str
            Binning strategy. Accepts two options: 'bins' and 'quantiles'. Defaults to 'bins'.
            'bins': input arrays are splitted into bins with equal
                and fixed steps based on 'expected' array
            'quantiles': input arrays are binned according to 'expected' array
                with given number of n_bins
            'beta': input arrays are uniformly binned within the interval of [0, 1]
                corresponding the possible methylation beta-values
        n_bins : int
            Number of buckets for binning. Defaults to 16.

    Returns
    -------
        A single float number
    """
    expected_freqs, observed_freqs = _compute_freqs(expected, observed, nbins, bucket_type)
    # Clip freaquencies to avoid zero division
    expected_freqs = np.clip(expected_freqs, a_min=a_min, a_max=None)
    observed_freqs = np.clip(observed_freqs, a_min=a_min, a_max=None)
    # Calculate PSI
    psi_value = (expected_freqs - observed_freqs) * np.log(expected_freqs / observed_freqs)
    # bound maximal psi terms
    psi_value = np.clip(psi_value, a_max=p_max, a_min=None)
    psi_value = np.sum(psi_value)
    return psi_value

def _jsd(expected, observed, nbins=16, bucket_type='bins'):    
    """Calculate Jensen-Shannon divergence metric for two arrays.   
    Parameters
    ----------
        expected : list-like
            Array of expected values
        actual : list-like
            Array of actual values
        bucket_type : str
            Binning strategy. Accepts two options: 'bins' and 'quantiles'. Defaults to 'bins'.
            'bins': input arrays are splitted into bins with equal
                and fixed steps based on 'expected' array
            'quantiles': input arrays are binned according to 'expected' array
                with given number of n_bins
            'beta': input arrays are uniformly binned within the interval of [0, 1]
                corresponding the possible methylation beta-values
        n_bins : int
            Number of buckets for binning. Defaults to 16.
    """
    from scipy.spatial.distance import jensenshannon
    expected_freqs, observed_freqs = _compute_freqs(expected, observed, nbins, bucket_type)
    return jensenshannon(expected_freqs, observed_freqs)

def _tvd(expected, observed, nbins=16, bucket_type='bins'):
    """Calculate Total Variation distance metric for two arrays.
    Parameters
    ----------
        expected : list-like
            Array of expected values
        actual : list-like
            Array of actual values
        bucket_type : str
            Binning strategy. Accepts two options: 'bins' and 'quantiles'. Defaults to 'bins'.
            'bins': input arrays are splitted into bins with equal
                and fixed steps based on 'expected' array
            'quantiles': input arrays are binned according to 'expected' array
                with given number of n_bins
            'beta': input arrays are uniformly binned within the interval of [0, 1]
                corresponding the possible methylation beta-values
        n_bins : int
            Number of buckets for binning. Defaults to 16.
    """
    expected_freqs, observed_freqs = _compute_freqs(expected, observed, nbins, bucket_type)
    return np.sum(np.abs(expected_freqs - observed_freqs))

def _ks(expected, observed):
    """Calculate Kolmogorov-Smirnov distance metric and p-value for two arrays.
       The zero p-value is replaced with a small constant.
    Parameters
    ----------
        expected : list-like
            Array of expected values
        actual : list-like
            Array of actual values
    """
    from scipy.stats import ks_2samp
    stat, pval = ks_2samp(expected, observed)
    pval = pval if pval > 1e-100 else 1e-100
    return stat, pval


def calculate_shift(
        expected: np.ndarray, 
        actual: np.ndarray, 
        metric_type: str = 'psi',
        bucket_type: str = "bins", 
        n_bins: int = 16,
) -> np.ndarray:
    """Apply PSI calculation to 2 1-d or 2-d arrays.

    Parameters
    ----------
    expected : list-like
        Array of expected values
    actual : list-like
        Array of actual values
    metric_type : str
        Type of metric measuring shift in a given covariate. Defaults to 'psi'.
            'psi' - population stability index 
            'jsd' - Jensen-Shannon divergence
            'tvd' - total variation distance
            'ks' - Kolmogorov-Smirnov statistic + p-value
    bucket_type : str
        Binning strategy. Accepts two options: 'bins' and 'quantiles'. Defaults to 'bins'.
            'bins' - input arrays are splitted into bins with equal
                and fixed steps based on ’expected' array
            'quantiles' - input arrays are binned according to ’expected’ array
                with given number of n_bins
    n_bins : int
        Number of buckets for binning. Defaults to 10.

    Returns
    -------
        np.ndarray of shift values
        in case of 'ks' also returns p-values
    """
    m = expected.shape[1]
    shift_values = np.empty(m)
    p_values = np.empty(m)
    for i in range(0, m):
        if metric_type == 'psi':
            shift_values[i] = _psi(expected[:, i], actual[:, i], n_bins, bucket_type)
        elif metric_type == 'jsd':
            shift_values[i] = _jsd(expected[:, i], actual[:, i], n_bins, bucket_type)
        elif metric_type == 'tvd':
            shift_values[i] = _tvd(expected[:, i], actual[:, i], n_bins, bucket_type)
        elif metric_type == 'ks':
            kstat, kpval = _ks(expected[:, i], actual[:, i])
            shift_values[i] = kstat
            p_values[i] = kpval
        else:
            raise NotImplementedError()
    if metric_type != 'ks':
        return shift_values
    else:
        return shift_values, p_values


def inverse_train_test_procedure(model, 
                                 X_train, y_train, 
                                 X_test, y_test=None, 
                                 params={},
                                 verbose=1):
    """
    Perform inversed train test procedure (ITTP) detecting the interchangeability between train and test datasets.
    The idea of the procedure is based on the reasonable assumption that once we obtained a good predictions for
    test datasets, the predictions themselves can be used for traning a new model. Which in its turn can properly
    predict the original train dataset.
    
    Parameters
    ----------
    model : sklearn.BaseModel
        A base model using for training on both steps. Can be any scikit-learn model.
        In the present research we used LassoCV because of its automatic hyperparameter search 
        and built-in CV procedure.
    
    X_train : pandas.DataFrame
        Train set - will be used for train on step 1, and for test on step 2.
    
    X_test : pandas.DataFrame
        Test set - will be used for test on step 1, and for train on step 2.
    
    y_train : pandas.Series or numpy.ndarray
        Train targets - will be used for train on step 1, and for test on step 2.

    y_test : pandas.Series or numpy.ndarray or None
        Test targets - if provided will be used for test on step 1, and for train on step 2 
        (falsifiable case of ITTP). If not provided, performance computation for test at 
        step 2 is impossible (unfalsifiable case of ITTP). 
    
    params : dict
        Any parameters you want to pass to the base model.
    

    Returns
    -------
        (nd.array, nd.array, nd.array, nd.array)
        first two are predictions on step1 train and test;
        last two are predictions on step2 train and test.
    """
    from sklearn.metrics import r2_score, mean_absolute_error
    # For Aging train/test sets
    #step 1: usual model fit
    model1 = model(**params)
    model1.fit(X_train, y_train)

    y_train_predict_step1 = model1.predict(X_train)
    y_test_predict_step1 = model1.predict(X_test)

    print('Step 1 results:')
    r2_tr_step1 = round(r2_score(y_train, y_train_predict_step1), 3)
    mae_tr_step1 = round(mean_absolute_error(y_train, y_train_predict_step1), 3)
    print(f'R2 train = {r2_tr_step1}')
    print(f'MAE train = {mae_tr_step1}')
    if y_test is not None:
        r2_te_step1 = round(r2_score(y_test, y_test_predict_step1), 3)
        mae_te_step1 = round(mean_absolute_error(y_test, y_test_predict_step1), 3)
        print(f'R2 test = {r2_te_step1}')
        print(f'MAE test = {mae_te_step1}')
    print('-'*30)

    #step 2: use the predicted y as response for fitting model on test, then predict train
    model2 = model(**params)
    model2.fit(X_test, y_test_predict_step1)

    y_train_predict_step2 = model2.predict(X_test)
    y_test_predict_step2 = model2.predict(X_train)

    print('Step 2 results:')
    if y_test is not None:
        r2_tr_step2 = round(r2_score(y_test, y_train_predict_step2), 3)
        mae_tr_step2 = round(mean_absolute_error(y_test, y_train_predict_step2), 3)
        print(f'R2 train = {r2_tr_step2}')
        print(f'MAE train = {mae_tr_step2}')
    
    r2_te_step2 = round(r2_score(y_train, y_test_predict_step2), 3)
    mae_te_step2 = round(mean_absolute_error(y_train, y_test_predict_step2), 3)
    print(f'R2 test = {r2_te_step2}')
    print(f'MAE test = {mae_te_step2}')
    print('-'*30) 
    
    return (y_train_predict_step1, 
            y_test_predict_step1, 
            y_train_predict_step2,
            y_test_predict_step2)