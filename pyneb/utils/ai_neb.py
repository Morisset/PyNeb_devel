# -*- coding: utf-8 -*-
"""
ANN regression manager used by PyNeb.

This is a simplified version of the original ai4neb manage_RM class,
restricted to scikit-learn MLPRegressor models.


@author: ChristopheMorisset & RogelioOrozcoDuarte
"""

# coding: utf-8
import numpy as np
import time
import random
from glob import glob

from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPRegressor



try:
    import joblib
    
except:
    try:
        from sklearn.externals import joblib
    except:
        pass





class manage_RM(object):
    """
    Regression model based on scikit-learn MLPRegressor
    """

    def __init__(self,
                 X_train=None, y_train=None, 
                 X_test=None, y_test=None,
                 scaling=False, 
                 use_log=False,
                 verbose=False,
                 RM_filename=None, 
                 random_seed=None,
                **kwargs):
        """
        Object to manage a Regression Model based on a scikit-learn ANN (MLPRegressor).
        
        Parameters
        ----------
        X_train, y_train : array-like
            Training datasets, input and output respectively.
        
        X_test, y_test : array-like
            Test datasets used for predictions.
        
        scaling : bool
            If True, input data are scaled using StandardScaler.
        
        use_log : bool
            If True, log10 is applied to the input data before scaling.
        
        RM_filename : str
            Filename used to load a previously saved model.
        
        random_seed : int
            Random seed used to initialize the ANN.
        """

        self.verbose = verbose
        self.init_random(random_seed)
        self.scaling = scaling
        self.use_log = use_log
        self.scaler = None
        self.train_scaled = False
        self.test_scaled = False
        self.isfin = None
        #self.RM_type = "SK_ANN"
        self.X_train = self._copy_None(X_train)
        self.y_train = self._copy_None(y_train)
        self.X_test = self._copy_None(X_test)
        self.y_test = self._copy_None(y_test)
        self._init_dims(train=True, test=True)
        self.RM_version = "0.17"
        if self.verbose:
            #check N_train and N_test
            print('Training set size = {}, Test set size = {}'.format(self.N_train, self.N_test))

        self.y_train_ori = self.y_train
        self.y_test_ori = self.y_test            
        if self.scaling:
            self.scale_sets(use_log=self.use_log)
        else:
            self.X_train_unscaled = self.X_train
            self.X_test_unscaled = self.X_test
            self.y_train_unscaled = self.y_train
        self.RMs = None
        self.trained = False
        self._multi_predic = True
        self.RM_filename = RM_filename

        if self.RM_filename is not None:
            self.load_RM(filename=self.RM_filename)

        if self.verbose:
            print('Training set size = {}, Test set size = {}'.format(self.N_train, self.N_test))
        
    def init_random(self, seed):
        self.random_seed = seed
        np.random.seed(self.random_seed)
        random.seed(self.random_seed)
        
        
    def _init_dims(self, train=True, test=True):
        
        def get_shape(a):
            if a is None:
                return 0, 0
            else:
                return a.shape
        if train:
            self.N_train, self.N_in = get_shape(self.X_train) 
            self.N_train_y, self.N_out = get_shape(self.y_train)
        if test:    
            self.N_test, self.N_in_test = get_shape(self.X_test)
            self.N_test_y, self.N_out_test = get_shape(self.y_test)
        
        
    def _copy_None(self, v, add_dim=True):
        if v is None:
            to_return = None
        else:
            to_return = v.copy()
            if np.ndim(v) == 1 and add_dim:
                to_return = np.expand_dims(to_return, axis=1)
        return to_return
    
    def init_RM(self, user_defined_RM=None, **kwargs):
        """
        Initialisation of the Regression Model.
        user_defined_RM: an optional dictionnary containing a (list of) model(s),
            train_params, and _multi_predic.
        Any parameter is passed to the Model.
        self.N_out RM can be needed if not ANN type. 
        They are stored in self.RMs list.
        """
        

        if user_defined_RM is not None:
            if type(user_defined_RM) is not dict:
                raise TypeError('user_defined_RM needs to be a dictionnary')
            if 'model' not in user_defined_RM:
                raise ValueError('user_defined_RM dictionnary does not contain model')
            else:
                user_defined_model = user_defined_RM['model']
            if type(user_defined_model) is list:
                self.RMs = user_defined_model
            else:
                self.RMs = [user_defined_model]
            if 'train_params' in user_defined_RM:
                self.train_params = user_defined_RM['train_params']
            if '_multi_predic' in user_defined_RM:
                self._multi_predic = user_defined_RM['_multi_predic']
            return

        self.RMs = []
        self.train_params = {}
        self._multi_predic = False
        self.RMs = [MLPRegressor(random_state=self.random_seed, **kwargs)]
        self._multi_predic = True

        

    def set_test(self, X=None, y=None):

        self.X_test = self._copy_None(X)
        self.y_test = self._copy_None(y)
                
        self.test_scaled = False
        if self.scaling:
            self.scale_sets(use_log=self.use_log)
        else:
            self.X_test_unscaled = self.X_test
         
        self._init_dims(train=False, test=True)

    def _log_data(self, X, y):
        """
        Apply log10 to X. 
        Filter X and y to isfinite(X)
        Return filtered X and y
        """
        if X is None:
            return None, None
        else:
            n_keys = X.shape[1]
            with np.errstate(invalid='ignore', divide='ignore'):
                X = np.log10(X)
            self.isfin = np.isfinite(X).sum(1) == n_keys
            X = X[self.isfin]
            if y is not None:
                y = y[self.isfin]
            
            return X, y
        
    def _set_scaler(self, force=False):
        if (self.scaler is None) or force:
            self.scaler = StandardScaler()
            self.scaler.fit(self.X_train)


    def scale_sets(self, force=False, use_log=False):
        """
        A scaler is created in self.scaler if it does not already exist.
        (it may have been recovered from saved file)
        It is fit to self.X_train if it is just created. Otherwise is 
        is supposed to already have been fit. To restart fitting, set self.scaler to None.
        If use_log, the self.X_train and self.X_test data are transformed to log10
        and only finite data is kept.
        The scaler is applied to self.X_train and self.X_test if they exist.
        """
        self.y_train_unscaled = self.y_train
        log_str = ''
        pca_str = ''
        if (not self.train_scaled) or force:
            self.X_train_unscaled = self._copy_None(self.X_train)
            if use_log:
                self.X_train, self.y_train = self._log_data(self.X_train, self.y_train)
                log_str = 'Log10 applied. '
            if self.X_train is not None:
                self._set_scaler()
                self.X_train = self.scaler.transform(self.X_train)

                self._init_dims(train=True, test=False)
                self.train_scaled = True
            
            if self.verbose:
                print('Train data scaled. {}{}'.format(log_str,pca_str))
        
        self.y_test_unscaled = self.y_test
        log_str = ''
        pca_str = ''
        if (not self.test_scaled) or force:
            if self.X_test is not None:
                self.X_test_unscaled = self._copy_None(self.X_test)
                if use_log:
                    self.X_test, self.y_test = self._log_data(self.X_test, self.y_test)
                    log_str = 'Log10 applied. '
                if self.X_test is not None:
                    self.X_test = self.scaler.transform(self.X_test)
                    self._init_dims(train=False, test=True)
                    self.test_scaled = True
            
            if self.verbose:
                print('Test data scaled. {}{}'.format(log_str,pca_str))
                
        if self.verbose:
            print('Training set size = {}, Test set size = {}'.format(self.N_train, self.N_test))

    def train_RM(self, p_backend='loky', scoring=True):
        """
        Training the models.
        """
        start = time.time()
        if not self.train_scaled and self.verbose:
            print('WARNING: training data not scaled')
        self.train_score = []
        self.history = []
        if self.N_train != self.N_train_y:
            raise Exception('N_train {} != N_train_y {}'.format(self.N_train,
                            self.N_train_y))
        if self.verbose:
            print('Training {} inputs for {} outputs with {} data'.format(self.N_in,
                  self.N_out, self.N_train_y))
        if self._multi_predic:
            RM = self.RMs[0]
            if self.y_train.ndim == 1:
                y_train = np.ravel(self.y_train)
            elif self.y_train.ndim == 2 and self.y_train.shape[1] == 1:
                y_train = np.ravel(self.y_train)
            else:
                y_train = self.y_train

            history = RM.fit(self.X_train, y_train, **self.train_params)
            self.history = [history]
            
            if scoring:
                train_score = score(RM, self.X_train, y_train)
                self.train_score = [train_score]
                iter_str = '.'
                if self.verbose:
                    try:
                        iter_str = ', with {} iterations.'.format(RM.n_iter_)
                    except:
                        iter_str = '.'
                        #train score is not tha useful here, but it does not hurt to keep it.
                    print('RM trained{} Score = {:.3f}'.format(iter_str, train_score))
        else:
             pass

        self.trained = True
        end = time.time()
        self.training_time = end - start
        if self.verbose:
            for RM in self.RMs:
                print(RM)
            print('Training time {:.1f} s.'.format(self.training_time))


    
    def predict(self, scoring=False, reduce_by=None, **kwargs):
        """
        Compute the prediction using self.X_test
        Results are stored into self.pred
        if scoring, a score is computed comparing with self.y_test
        """
        start = time.time()
        if self.RMs is None:
            raise Exception('WARNING: Regression Model not set up')
        if not self.trained:
            raise Exception('WARNING: Regression Model not trained')
        if not self.test_scaled and self.verbose:
            print('WARNING: test data not scaled')
        if self._multi_predic:
            #this is the only line needed for multipredict
            self.pred = self.RMs[0].predict(self.X_test, **kwargs)

        self.pred = np.asarray(self.pred) # in case of a Tensor
        end = time.time()
        if self.verbose:
            print('Predicting from {} inputs to {} outputs using {} data in {:.2f} secs.'.format(self.N_in_test,
                  self.N_out, self.N_test, end - start))
        
    def save_RM(self, filename='RM', save_train=False, save_test=False, **kwargs):
        """
        Save the current regression model and its configuration.
        
        The following information is stored:
        
            RM_version, RM_type
            X_train, y_train, X_test, y_test
            scaling, use_log
            train_scaled, test_scaled
            scaler
            N_in, N_out, N_in_test, N_out_test
            N_test, N_test_y, N_train, N_train_y
            train_score
            _multi_predic
            trained, training_time
            random_seed
            RMs
        
        The model is saved using joblib in the file:
        
            filename.ai4neb_sk
        """
                
        if not self.trained:
            raise Exception('Regression Model not trained')
        if save_train:
            X_train, y_train = self.X_train_unscaled, self.y_train_unscaled
        else:
            X_train, y_train = None, None
        if save_test:
            X_test, y_test = self.X_test_unscaled, self.y_test_unscaled
        else:
            X_test, y_test = None, None
        

        to_save = [
                    self.RM_version,
                    X_train, y_train, X_test, y_test,
                    self.scaling, self.use_log,
                    self.train_scaled, self.test_scaled,
                    self.scaler,
                    self.N_in, self.N_out,
                    self.N_in_test, self.N_out_test,
                    self.N_test, self.N_test_y,
                    self.N_train, self.N_train_y,
                    self.train_score, self._multi_predic,
                    self.trained, self.training_time,
                    self.random_seed,
                    self.RMs
                ]

        joblib.dump(to_save, filename+".ai4neb_sk", **kwargs)
        # if self.RM_type[0:3] == 'SK_': 
        #     joblib.dump(to_save, filename+'.ai4neb_sk', **kwargs)
        #     if self.verbose:
        #         print('RM save to {}.ai4neb_sk'.format(filename))

        # else:
        #    print('Do not know how to save {} machine'.format(self.RM_type))
                    
    def load_RM(self, filename='RM', notry=False, compile_=False):
        """
        Loading previously saved model.
        joblib is used to load.
        A WARNING is issued if a different version is found in the file.
        
        Example: Only X_test is needed.
        
        RM = manage_RM(X_test=X_test)
        RM.load_RM(filename)
        
        it can also be included to the instantiation:
            
        RM = manage_RM(X_test=X_test, RM_filename=RM_filename)
        RM.scale_sets(use_log=True)
        RM.predict(scoring=False)
        
        """
        files = glob("{}.*".format(filename))
        
        
        if "{}.ai4neb_sk".format(filename) in files:
            to_read = "{}.ai4neb_sk".format(filename)
            

        else:
            to_read = None
            print('No ai4neb file found for {}'.format(filename))
            self.model_read = False
            return
        
        if notry:
            RM_tuple = joblib.load(to_read)
            if self.verbose:
                print('RM loaded from {}'.format(to_read))
        else:
            try:
                RM_tuple = joblib.load(to_read)
                if self.verbose:
                    print('RM loaded from {}'.format(to_read))
            except:
                print('!! ERROR reading {}'.format(to_read))
        
        load_version = RM_tuple[0]
        if self.RM_version != load_version and self.verbose:
            print('WARNING: version loaded from {} is {}. Version from RM class is {}.'.format(to_read, 
                                                                  load_version, self.RM_version))
        
        if load_version == "0.17":
            (self.RM_version, 
                   self.X_train, self.y_train, self.X_test, self.y_test,
                   self.scaling,  
                   self.use_log, 
                   self.train_scaled, self.test_scaled,
                   self.scaler,  
                   self.N_in, self.N_out, self.N_in_test, self.N_out_test,
                   self.N_test, self.N_test_y, self.N_train, self.N_train_y,
                   self.train_score, self._multi_predic,
                   self.trained, self.training_time,
                   self.random_seed, 
                   self.RMs) = RM_tuple
        else:
            self.model_read = False
            print('!! ERROR. This version is not supported.')
          
        if self.scaling:
            self.scale_sets(use_log=self.use_log, force=True)
        else:
            self.X_train_unscaled = self.X_train
            self.X_test_unscaled = self.X_test
            self.y_train_unscaled = self.y_train
        self.model_read =True
        
def score(RM, X, y_true, axis=None):
    """
    (1 - u/v), where u is the residual sum of squares ((y_true - y_pred) ** 2).sum() 
    and v is the total sum of squares ((y_true - y_true.mean()) ** 2).sum().
    """

    y_pred = RM.predict(X)
    if y_pred.ndim == 2 and y_pred.shape[1] == 1:
            y_pred = np.ravel(y_pred)
    u = ((y_true - y_pred) ** 2).sum(axis=axis)
    v = ((y_true - y_true.mean()) ** 2).sum(axis=axis)
    
    return 1 - u/v
                
#%% __main__
if __name__ == "__main__":
    pass