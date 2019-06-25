v1.4.1

This is a minor release

API Changes
~~~~~~~~~~~


Enhancements
~~~~~~~~~~~~
* Sky diffuse from `pvl_perez.m` is limited between 0.0 W/m^2 and extraterrestrial horizontal irradiance. Values outside this range are set to NaN.

Bug fixes
~~~~~~~~~
* Adjusted convergence criteria in `pvl_lambertw.m`. `pvl_singlediode.m` should converge more consistently and faster.


Contributors
~~~~~~~~~~~~
* Jonathan Allen (:ghuser:'joallen6')
