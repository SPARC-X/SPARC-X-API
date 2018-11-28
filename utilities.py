#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 12:52:34 2018

@author: benjamin
"""

import numpy as np
from ase.data import chemical_symbols
from scipy.misc import factorial


valences = [0]+[1,2]+[1,2,3,4,5,6,7,8]*2+ \
           [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]*2+ \
           list(np.arange(1,33))*2
           #atomic numbers correspond to traditional number of valence electrons
valence_dict = {}
for number,element in enumerate(chemical_symbols):
    valence_dict[element] = valences[number]

def h2gpts(h, cell_cv, idiv=4):
    """Convert grid spacing to number of grid points divisible by idiv.
    Taken from GPAW:
        https://gitlab.com/gpaw/gpaw/blob/master/gpaw/utilities/__init__.py

    Note that units of h and cell_cv must match!

    h: float
        Desired grid spacing in.
    cell_cv: 3x3 ndarray
        Unit cell.
    """

    L_c = (np.linalg.inv(cell_cv)**2).sum(0)**-0.5
    return np.maximum(idiv, (L_c / h / idiv + 0.5).astype(int) * idiv)

def Ecut2h(FDn,Ecut,tol=0.1):
    #to do complete python conversion and integrate
    """
    converted from MATLAB code provided by Qimen
    """
    w2 = np.zeros(FDn+1)
    for i in range(FDn+1):
        k = i+1
        w2[i+1] = (2*(-1)**(k+1))*(factorial(FDn)**2)/ \
                    (k*k*factorial(FDn-k)*factorial(FDn+k))
        w2[0] = w2[0]-2*(1/(k*k))
    kk = np.linspace(0,np.pi,1001)
    y_cos = -w2[0]+ (-2*w2[1:]) * np.cos(np.arange(1,FDn) * kk)
    """
    % Find correlation bw/ Ecut (Ry) and mesh size h (Bohr)
    clear all; close all; clc;
    
    % input variables
    FDn = 6; % FD order / 2
    Ecut = 30; % in Ry
    
    
    
    
    
    
    epsilon = 1e-1; % tolerance threshold for FD 2nd derivative approx.
    % Ecut in Ha
    %Ecut = Ecut * 0.5;
    
    % finite difference weights
    w2 = zeros(1,FDn+1); 
    for k=1:FDn
        w2(k+1) = (2*(-1)^(k+1))*(factorial(FDn)^2)/...
                        (k*k*factorial(FDn-k)*factorial(FDn+k));
        w2(1) = w2(1)-2*(1/(k*k));
    end
    
    kk = linspace(0,pi,1001);
    y_cos =  -w2(1) + (-2*w2(2:end)) * cos((1:FDn)' * kk);
    freq_err = abs(y_cos - kk.^2);
    kc = kk(max(find(freq_err < epsilon)))
    
    plot(kk,y_cos,'-',kk,kk.^2,'--','LineWidth',2);
    hold on
    plot(kc,-w2(1) + (-2*w2(2:end)) * cos((1:FDn)' * kc),'p','MarkerSize',15);
    hold off
    
    kc_by_pi = kc / pi;
    
    h = kc / sqrt(Ecut);
    
    fprintf('******************************************\n');
    fprintf('Ecut = %.2f Ry = %.2f Ha corresponds to:\n\t h = %.3f Bohr in real space\n', ...
            Ecut,0.5*Ecut,h);
    fprintf('******************************************\n');
    
    
    h_pw = pi / sqrt(Ecut)
    
    N_r = (h_pw / h)^3
    """
    
def h2Ecut():
    #to do, convert to python and test
    """
    % Find correlation bw/ Ecut (Ry) and mesh size h (Bohr)
    clear all; close all; clc;
    
    % input variables
    FDn = 6; % FD order / 2
    Ecut = 30; % in Ry
    
    
    
    
    
    
    epsilon = 1e-1; % tolerance threshold for FD 2nd derivative approx.
    % Ecut in Ha
    %Ecut = Ecut * 0.5;
    
    % finite difference weights
    w2 = zeros(1,FDn+1); 
    for k=1:FDn
        w2(k+1) = (2*(-1)^(k+1))*(factorial(FDn)^2)/...
                        (k*k*factorial(FDn-k)*factorial(FDn+k));
        w2(1) = w2(1)-2*(1/(k*k));
    end
    
    kk = linspace(0,pi,1001);
    y_cos =  -w2(1) + (-2*w2(2:end)) * cos((1:FDn)' * kk);
    freq_err = abs(y_cos - kk.^2);
    kc = kk(max(find(freq_err < epsilon)))
    
    plot(kk,y_cos,'-',kk,kk.^2,'--','LineWidth',2);
    hold on
    plot(kc,-w2(1) + (-2*w2(2:end)) * cos((1:FDn)' * kc),'p','MarkerSize',15);
    hold off
    
    kc_by_pi = kc / pi;
    
    h = kc / sqrt(Ecut);
    
    fprintf('******************************************\n');
    fprintf('Ecut = %.2f Ry = %.2f Ha corresponds to:\n\t h = %.3f Bohr in real space\n', ...
            Ecut,0.5*Ecut,h);
    fprintf('******************************************\n');
    
    
    h_pw = pi / sqrt(Ecut)
    
    N_r = (h_pw / h)^3
    """    
        
