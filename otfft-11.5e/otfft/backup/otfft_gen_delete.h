switch (log_N) {
case  0: break;
case  1: delete static_cast<FFT3*>(obj); break;
case  2: delete static_cast<FFT1*>(obj); break;
case  3: delete static_cast<FFT7*>(obj); break;
case  4: delete static_cast<FFT2*>(obj); break;
case  5: delete static_cast<FFT7*>(obj); break;
case  6: delete static_cast<FFT1*>(obj); break;
case  7: delete static_cast<FFT1*>(obj); break;
case  8: delete static_cast<FFT1*>(obj); break;
case  9: delete static_cast<FFT3*>(obj); break;
case 10: delete static_cast<FFT1*>(obj); break;
case 11: delete static_cast<FFT5*>(obj); break;
case 12: delete static_cast<FFT6*>(obj); break;
case 13: delete static_cast<FFT7*>(obj); break;
case 14: delete static_cast<FFT5*>(obj); break;
case 15: delete static_cast<FFT6*>(obj); break;
case 16: delete static_cast<FFT6*>(obj); break;
case 17: delete static_cast<FFT4*>(obj); break;
case 18: delete static_cast<FFT7*>(obj); break;
case 19: delete static_cast<FFT7*>(obj); break;
case 20: delete static_cast<FFT5*>(obj); break;
case 21: delete static_cast<FFT7*>(obj); break;
case 22: delete static_cast<FFT6*>(obj); break;
case 23: delete static_cast<FFT6*>(obj); break;
case 24: delete static_cast<FFT5*>(obj); break;
default: delete static_cast<FFT8*>(obj); break;
}
