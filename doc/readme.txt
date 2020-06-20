License for research purpose only. Expires on 2020-12-31. Contact Xiaoming Chen (chenxiaoming@ict.ac.cn) for a new license or download a new version once expired.

Three demos are provided. demo_real.c and demo_complex.c demonstrate basic usage of NICSLU for real and complex numbers respectively. demo.c demonstrates a more practical usage of NICSLU. 
For Windows users: to run the demos, please put "lib32/nicslu.dll" and "license/nicslu.lic" together with the executable files.
For Linux users: to run the demos, please put "license/nicslu.lic" together with the executable files.

Please read faq.pdf for some very common questions. These questions help you understand how to use NICSLU and some potential problems you may meet.

To use lib64_int64, you must define the macro __NICS_INT64. Please define the macro __NICS_INT64 by the compiler (this will ensure that this macro is seen by any code), instead of in any C code.

NICSLU DLLs require Windows 7 or higher version Windows.
