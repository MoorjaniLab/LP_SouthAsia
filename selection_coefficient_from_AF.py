#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  5 13:44:03 2025

@author: elisekerdoncuff
"""

import math

def solve_s_additive(p0, pt, t):
    return (2/t) * math.log((pt*(1-p0)) / (p0*(1-pt)))

def solve_s_dominant(p0, pt, t):
    return (1/t) * (math.log((pt*(1-p0)) / (p0*(1-pt))) + 1/(1-pt) -1/(1-p0) )

# Toda example usage:
p0 = 0.078       # initial allele freq
p_obs = 0.6111   # observed after dt generations
dt = 125          # number of generations between samples




s_dominant = solve_s_dominant(p0, p_obs, dt, s_min=-1, s_max=1)
s_additive = solve_s_additive(p0, p_obs, dt, s_min=-1, s_max=1)
print("Estimated s_dominant =", s_dominant)
print("Estimated s_additive =", s_additive)
