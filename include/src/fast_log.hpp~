#ifndef FAST_LOG_HPP
#define FAST_LOG_HPP

//fast log is copied from http://www.flipcode.com/archives/Fast_log_Function.shtml
//to see how floats are stored in memory see http://softwareengineering.stackexchange.com/questions/215065/can-anyone-explain-representation-of-float-in-memory
inline float fast_log2 (float val){
  int * const    exp_ptr = reinterpret_cast<int*>(&val);
  int            x = *exp_ptr;  //this is casting val to integer type
  const int      log_2 = ((x >> 23) & 255) - 128; //this is calculating the "exponent - 1", which is stored in bits 2-9 of the float (first bit is sign)
  x &= ~(255 << 23);
  x += 127 << 23; //x is now the mantissa which is a number between 1 and 2
  *exp_ptr = x; //reinterpret x as float

  val = ((-1.0f/3) * val + 2) * val - 2.0f/3;   // (1) //computes 1+log2(val) using a polynomial approximation

  return (val + log_2);
}

inline float fast_log (const float &val){
  return (fast_log2 (val) * 0.69314718f);
} 


#endif //DATA_HPP
