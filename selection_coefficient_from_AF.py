#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  5 13:44:03 2025

@author: elisekerdoncuff
"""


import numpy as np
from scipy.optimize import brentq
import math

def f_onegen_dominant(p, s):
    # one-generation deterministic update for h=1 (dominant)
    return (p*(1+s)) / (1 + s*(2*p - p*p))

def iterate_map_dominant(p0, s, dt):
    p = p0
    for _ in range(dt):
        p = f_onegen_dominant(p, s)
    return p

def solve_s_from_p0_dominant(p0, p_obs, dt, s_min=-100.0, s_max=100.0):
    # root function: for a given s, iterated map minus observed frequency
    def rootfun(s):
        return iterate_map_dominant(p0, s, dt) - p_obs
    # bracket and solve
    return brentq(rootfun, s_min, s_max)


def iterate_additive(p0, s, t):
    p = p0
    for _ in range(t):
        wAA = 1+s
        wAa = 1+s/2
        waa = 1
        mean_w = p**2*wAA + 2*p*(1-p)*wAa + (1-p)**2*waa
        p = (p**2*wAA + p*(1-p)*wAa)/mean_w
    return p

def solve_s_from_p0_additive(p0, pt, t, s_min=0, s_max=1):
    from scipy.optimize import brentq
    def rootfun(s):
        return iterate_additive(p0, s, t) - pt
    return brentq(rootfun, s_min, s_max)

# Toda example usage:
p0 = 0.078       # initial allele freq
p_obs = 0.6111   # observed after dt generations
dt = 125          # number of generations between samples




s_dominant = solve_s_from_p0_dominant(p0, p_obs, dt, s_min=-1, s_max=1)
s_additive = solve_s_from_p0_additive(p0, p_obs, dt, s_min=-1, s_max=1)
print("Estimated s_dominant =", s_dominant)
print("Estimated s_additive =", s_additive)