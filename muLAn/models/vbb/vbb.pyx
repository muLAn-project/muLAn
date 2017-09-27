# -*- coding: utf-8 -*-

cdef extern from "VBBinaryLensingLibrary.cpp":
    double magVBBU(double s, double q, double rho, double y1, double y2, double accuracy)
    double magVBBLLD(double s, double q, double rho, double gamma, double y1, double y2, double accuracy)

def vbbmagU(s, q, rho, y1, y2, accuracy):
    return magVBBU(s, q, rho, y1, y2, accuracy)

def vbbmagLLD(s, q, rho, gamma, y1, y2, accuracy):
    return magVBBLLD(s, q, rho, gamma, y1, y2, accuracy)

