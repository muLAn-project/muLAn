# -*- coding: utf-8 -*-

cdef extern from "VBBinaryLensingLibrary.cpp":
    double magVBB(double s, double q, double rho, double gamma, double y1, double y2, double accuracy)

def contourmag(s, q, rho, gamma, y1, y2, accuracy):
    return magVBB(s, q, rho, gamma, y1, y2, accuracy)

