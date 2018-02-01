#-------------------------------------------------------------------------------
# Name:        Generalized Logistic Regression
# Purpose:
#
# Author:      Michael Peter McGurk
#
# Created:     10/01/2017
# Copyright:   (c) Michael Peter McGurk 2017
#-------------------------------------------------------------------------------

import numpy
import scipy

import sklearn
import sklearn.linear_model
from sklearn import gaussian_process
from sklearn.base import BaseEstimator

class GeneralizedLogisticRegressor(BaseEstimator):
    def __init__(self, K=1.):
        self.K=K
        pass

    def fit(self,X,y):
        self.A=min(y)
        min_y=numpy.min(y)
        max_y=numpy.max(y)
        x_array=numpy.array(X)
        y_array=numpy.array(y)
        M_init=numpy.median( x_array[(y_array!=min_y)*(y_array!=max_y)])
        p0=[1,1,M_init,1]
        print p0
        self.best_param, self.cov=scipy.optimize.curve_fit(self.GenLog, X,y, p0=p0)
        self.B, self.Q, self.M, self.v=self.best_param

    def predict(self,X,y=None):
        return self.GenLog(X, self.B, self.Q, self.M, self.v)

    def score(self, X,y):
        predict=self.predict(X)
        return numpy.mean( (y-predict)**2)**.5

    def GenLog(self,x, B, Q, M, v):
        """Generalized logistic function."""
        C=1
        K=1.
        alpha=B*(x-M)
        return self.A+ (self.K-self.A)/(C+Q*numpy.exp(-1*alpha))**(1/v)

class JanoschekRegressor(BaseEstimator):
    def __init__(self, A=.5, K=1.):
        self.upper=K
        self.K=K
        pass

    def fit(self,X,y):
        self.best_param, self.cov=scipy.optimize.curve_fit(self.Janoschek, X,y,p0=[0,1,1])
        self.lower, self.growth_rate, self.inflection_point=self.best_param

    def predict(self,X,y=None):
        return self.Janoschek(X, self.lower, self.growth_rate, self.inflection_point)

    def score(self, X,y):
        predict=self.predict(X)
        return numpy.mean( (y-predict)**2)**.5
    def Janoschek(self,x,l,g, i):
        """Generalized logistic function."""


        return self.upper-(self.upper-l)*numpy.exp(-g*x**i)


class DevGeneralizedLogisticRegressor(BaseEstimator):
    def __init__(self, A=.5, K=1.):
        self.A=A
        self.K=K
        pass

    def fit(self,X,y):
        init=X[ numpy.argmin(y)]
        self.best_param, self.cov=scipy.optimize.curve_fit(self.DevGenLog, X,y,p0=[1,1e-3,init,1])
        self.C, self.B, self.M, self.v=self.best_param

    def predict(self,X,y=None):
        return self.DevGenLog(X, self.C, self.B, self.M, self.v)

    def score(self, X,y):
        predict=self.predict(X)
        return numpy.mean( (y-predict)**2)**.5

    def DevGenLog(self,x,C ,B, M, v):
        """Generalized logistic function."""

        alpha=numpy.exp(-1* B*(x-M))
        return 1-(C/v)*alpha*(alpha+1)**(-1./v-1)


class FMRegressor(BaseEstimator):
    def __init__(self, A=.5, K=1.):
        self.A=A
        self.K=K
        pass

    def fit(self,X,y):
        init=X[ numpy.argmin(y)]
        self.best_param, self.cov=scipy.optimize.curve_fit(self.FMReg, X,y,p0=[1,1,1,1,1, 1,1e-2,init,1])
        self.A1, self.B1,self.Q1,self.M1,self.v1, self.C2,self.B2,self.M2,self.v2=self.best_param

    def predict(self,X,y=None):
        return self.FMReg(X,self.A1, self.B1,self.Q1,self.M1,self.v1, self.C2,self.B2,self.M2,self.v2)

    def score(self, X,y):
        predict=self.predict(X)
        return numpy.mean( (y-predict)**2)**.5

    def DevGenLog(self,x,C ,B, M, v):
        """Generalized logistic function."""

        alpha=numpy.exp(-1* B*(x-M))
        return 1-(C/v)*alpha*(alpha+1)**(-1./v-1)

    def GenLog(self,x,A, B, Q, M, v):
        """Generalized logistic function."""
        C=1
        K=1.
        alpha=B*(x-M)
        return self.A+ (self.K-self.A)/(C+Q*numpy.exp(-1*alpha))**(1/v)

    def FMReg(x,A1, B1,Q1,M1,v1, C2,B2,M2,v2):
        Prec=GenLog(x,A1,B1,Q1,M1,v1)
        Rec=DevGenLog(x,C2,B2,M2,v2)
        return (Prec*Rec)**.5

def main():
    pass

if __name__ == '__main__':
    main()
