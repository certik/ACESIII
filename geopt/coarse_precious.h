/*
 * This should list the minimum number of files needed to restart a geometry
 * optimization or finite difference calculation.
 */
const char *coarse_precious[] = {
 "ZMAT.BAS", /* basis set cache */
 "JOBARC", /* job archive */
 "JAINDX", /* JOBARC metadata */
 "OPTARC", /* geom opt archive */
 "OPTARCBK", /* true geom opt archive in numerical geom opt calc */
 "DIPDER", /* dipole derivatives for vib freq intensities */
 "" /* must be NULL terminated */
};
