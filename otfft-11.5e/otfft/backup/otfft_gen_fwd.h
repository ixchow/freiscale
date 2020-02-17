switch (log_N) {
case  0: break;
case  1: static_cast<FFT3*>(obj)->fwd(x, y); break;
case  2: static_cast<FFT1*>(obj)->fwd(x, y); break;
case  3: static_cast<FFT7*>(obj)->fwd(x, y); break;
case  4: static_cast<FFT2*>(obj)->fwd(x, y); break;
case  5: static_cast<FFT7*>(obj)->fwd(x, y); break;
case  6: static_cast<FFT1*>(obj)->fwd(x, y); break;
case  7: static_cast<FFT1*>(obj)->fwd(x, y); break;
case  8: static_cast<FFT1*>(obj)->fwd(x, y); break;
case  9: static_cast<FFT3*>(obj)->fwd(x, y); break;
case 10: static_cast<FFT1*>(obj)->fwd(x, y); break;
case 11: static_cast<FFT5*>(obj)->fwd(x, y); break;
case 12: static_cast<FFT6*>(obj)->fwd(x, y); break;
case 13: static_cast<FFT7*>(obj)->fwd(x, y); break;
case 14: static_cast<FFT5*>(obj)->fwd(x, y); break;
case 15: static_cast<FFT6*>(obj)->fwd(x, y); break;
case 16: static_cast<FFT6*>(obj)->fwd(x, y); break;
case 17: static_cast<FFT4*>(obj)->fwd(x, y); break;
case 18: static_cast<FFT7*>(obj)->fwd(x, y); break;
case 19: static_cast<FFT7*>(obj)->fwd(x, y); break;
case 20: static_cast<FFT5*>(obj)->fwd(x, y); break;
case 21: static_cast<FFT7*>(obj)->fwd(x, y); break;
case 22: static_cast<FFT6*>(obj)->fwd(x, y); break;
case 23: static_cast<FFT6*>(obj)->fwd(x, y); break;
case 24: static_cast<FFT5*>(obj)->fwd(x, y); break;
default: static_cast<FFT8*>(obj)->fwd(x, y); break;
}
