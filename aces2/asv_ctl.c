/*
*  Copyright (c) 2003-2010 University of Florida
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  The GNU General Public License is included in this distribution
*  in the file COPYRIGHT.
*/ 
#include <unistd.h>	/* for NULL */
#include <string.h>	/* for strcmp */
#include <strings.h>	/* for strcasecmp */
#include <ctype.h>	/* for isspace */
#include <stdlib.h>	/* for atof, atol, getenv */
#include <stdio.h>	/* for fflush */

#ifndef NULL
#define NULL 0
#endif
#include "f77_name.h"	/* for F77_NAME */
#include "f_types.h"	/* for f_int */

extern void F77_NAME(errex,ERREX)();

/******************************************************************************/

/* a private value that the program can switch off with asv_hush() */
short bVerbose_Update = 1;

/* should we print tons of junk during the ASV updates? */
/*
#define PRINT_ALL
*/

#ifdef PRINT_ALL
void F77_NAME(asv_hush,ASV_HUSH)()
{ return; }
#else
void F77_NAME(asv_hush,ASV_HUSH)()
{ bVerbose_Update=0; return; }
#endif /* PRINT_ALL */

/******************************************************************************/

/* define type handles */
#define h_ICHAR_handle		0
#define h_ICHAR_string		1
#define h_ICHAR_f_int		2
#define h_ICHAR_double		3
#define h_ICHAR_f_int_array	4
#define h_ICHAR_double_array	5
const char *ichar_types[] =
{
    "handle",
    "string",
    "f_int",
    "double",
    "f_int array",
    "double array"
};

/******************************************************************************/

#define MONSTER_FLAGS /* use ioppar[600] instead of iflags[100], iflags2[500] */
#include "jodaflags.com" /* for f_flags.ioppar */

/* define the ASV namelist storage structure */
typedef struct ASV_nl_t
{
    char *alias;  /* a plain string for whole key token matching */
    char *oppar;  /* a key definition with a stub delimiter */
    f_int ichar;  /* a value type */
    f_int ideflt; /* the default value */
    char *units;  /* value units (for printing) */
} ASV_nl_t;

/*
 * ---=== RULES FOR ADDING A NEW ACES STATE VARIABLE ===---
 *
 * 1) DO NOT rearrange the order of the ASV_nl array unless you are prepared
 *    to edit every single line of source in ACES that requires an iflags
 *    or iflags2 lookup. The order of the definitions is critical to
 *    addressing them.
 *
 *    a) If an ASV is obsolete and must be removed, then blank out the whole
 *       line with NULL strings and zeroes.
 *
 *    b) If an ASV must be changed, run a test calculation with all possible
 *       ways to convert the key token. This will ensure the change does
 *       not conflict with a current ASV. However, this will not ensure a
 *       current ASV will not conflict with your change. BE CAREFUL!
 *
 * 2) Regardless of the type of value associated with your ASV, be sure
 *    to add your index to the defines in flags.h.
 *
 * 3) If your ASV accepts value tokens of type 'handle', then you must
 *    edit the asv_handle_proc routine. Instructions can be found there.
 *
 * 4) If your ASV accepts value tokens of type 'string' or any array, then
 *    make sure the largest (sized) value is smaller than the size of the
 *    work array that is passed in. For example, if your ASV accepts
 *    a f_int array of dimensions 8x3 (irreps,spin_pairs) then the data
 *    array must be at least 24 (Fortran) integers long.
 *
 */

/* define the ASVs */
#define MAX_KEY_LENGTH 16 /* including stub delimiter and NULL */
#define STUB_DELIM "#"
const ASV_nl_t ASV_nl[] =
{
/* COMMENTED NUMBERS ARE FORTRAN INDICES !!! */
/*       *alias,	*oppar,		ichar,		ideflt,	*units */
{/*1*/   "IPRNT",	"PRINT#",	h_ICHAR_f_int,	0,	""},
{/*2*/   "CALTYPE",	"CALC#LEVEL",	h_ICHAR_handle,	0,	""},
{/*3*/   "IDRLVL",	"DERIV#_LEV",	h_ICHAR_handle,	-1,	""},
{/*4*/   "ICCCNV",	"CC_C#ONV",	h_ICHAR_f_int,	7,	"(tol)"},
{/*5*/   "ISCFCV",	"SCF_C#ONV",	h_ICHAR_f_int,	7,	"(tol)"},
{/*6*/   "IXFTOL",	"XFORM#_TOL",	h_ICHAR_f_int,	11,	"(tol)"},
{/*7*/   "ICCCYC",	"CC_MAX#CYC",	h_ICHAR_f_int,	0,	"cycles"},
{/*8*/   "ILINDP",	"LIN#DEP_TOL",	h_ICHAR_f_int,	8,	""},
{/*9*/   "IRDOFM",	"RDO#",		h_ICHAR_handle,	-1,	""},
{/*10*/  "IRPP",	"SCF_EXTRAP#",	h_ICHAR_handle,	1,	""},
{/*11*/  "IREFNC",	"REF#ERENCE",	h_ICHAR_handle,	0,	""},
{/*12*/  "ICCEOR",	"CC_EXP#ORDER",	h_ICHAR_f_int,	0,	""},
{/*13*/  "IEVERY",	"TAM#P_SUM",	h_ICHAR_f_int,	0,	""},
{/*14*/  "ITOPT2",	"NTO#P_TAMP",	h_ICHAR_f_int,	15,	""},
{/*15*/  "ISCFDP",	"DAMPSCF#",	h_ICHAR_f_int,	20,	"x 0.01"},
{/*16*/  "ISCFCY",	"SCF_M#AXCYC",	h_ICHAR_f_int,	150,	"cycles"},
{/*17*/  "IOCCU",	"OCC#UPATION",	h_ICHAR_string,	0,	""},
{/*18*/  "IPROPS",	"PROP#S",	h_ICHAR_handle,	0,	""},
{/*19*/  "IDENS",	"DENS#ITY",	h_ICHAR_handle,	0,	""},
{/*20*/  "IRPPOR",	"SCF_EXPOR#DE",	h_ICHAR_f_int,	6,	""},
/*       *alias,	*oppar,		ichar,		ideflt,	*units */
{/*21*/  "ICCEXT",	"CC_EXT#RAPOL",	h_ICHAR_handle,	1,	""},
{/*22*/  "IBRKNR",	"BRUEC#KNER",	h_ICHAR_handle,	0,	""},
{/*23*/  "IXEFLD",	"XFI#ELD",	h_ICHAR_f_int,	0,	"x 10-6"},
{/*24*/  "IYEFLD",	"YFI#ELD",	h_ICHAR_f_int,	0,	"x 10-6"},
{/*25*/  "IZEFLD",	"ZFI#ELD",	h_ICHAR_f_int,	0,	"x 10-6"},
{/*26*/  "ISVINT",	"SAVE#_INTS",	h_ICHAR_handle,	0,	""},
{/*27*/  "IDRPMO",	"DROP#MO",	h_ICHAR_string,	0,	""},
{/*28*/  "ICHRGE",	"CHARGE#",	h_ICHAR_f_int,	0,	""},
{/*29*/  "IMULTP",	"MULT#IPLICTY",	h_ICHAR_f_int,	1,	""},
{/*30*/  "ICPHFT",	"CPHF_C#ONVER",	h_ICHAR_f_int,	12,	"(tol)"},
{/*31*/  "ICPHFC",	"CPHF_M#AXCYC",	h_ICHAR_f_int,	64,	"cycles"},
{/*32*/  "","",0,0,""},
{/*33*/  "","",0,0,""},
{/*34*/  "IQRHFO",	"QRHF_O#RBITA",	h_ICHAR_string,	0,	"offset"},
{/*35*/  "INCORE",	"INC#ORE",	h_ICHAR_handle,	0,	""},
{/*36*/  "IMEMSZ",	"MEM#ORY_SIZE",	h_ICHAR_f_int,	15000000,"Words"},
{/*37*/  "IFLREC",	"FILE_REC#SIZ",	h_ICHAR_f_int,	-1,	"Words"},
{/*38*/  "NON-HF",	"NONHF#",	h_ICHAR_handle,	0,	""},
{/*39*/  "IORBTP",	"ORB#ITALS",	h_ICHAR_handle,	-1,	""},
{/*40*/  "IRPPLS",	"SCF_EXPST#AR",	h_ICHAR_f_int,	8,	""},
/*       *alias,	*oppar,		ichar,		ideflt,	*units */
{/*41*/  "ILOCOC",	"LOCK_ORBO#CC",	h_ICHAR_handle,	0,	""},
{/*42*/  "ISTRIP",	"FILE_STRI#PE",	h_ICHAR_f_int,	0,	""},
{/*43*/  "HBAR",	"DOHBAR#",	h_ICHAR_handle,	0,	""},
{/*44*/  "ICHREC",	"CACHE#_RECS",	h_ICHAR_f_int,	-1,	""},
{/*45*/  "IGUESS",	"GUE#SS",	h_ICHAR_handle,	0,	""},
{/*46*/  "IJPRNT",	"JOD#A_PRINT",	h_ICHAR_f_int,	0,	""},
{/*47*/  "INR",		"OPT_METHO#D",	h_ICHAR_handle,	0,	""},
{/*48*/  "ICONTL",	"CONV#ERGENCE",	h_ICHAR_f_int,	4,	"H/bohr"},
{/*49*/  "IVEC",	"EIG#ENVECTOR",	h_ICHAR_f_int,	1,	""},
{/*50*/  "IDIE",	"NEG#EVAL",	h_ICHAR_handle,	2,	""},
{/*51*/  "ICURVY",	"CUR#VILINEAR",	h_ICHAR_handle,	0,	""},
{/*52*/  "ISTCRT",	"STP_SIZ_CTL#",	h_ICHAR_handle,	0,	""},
{/*53*/  "IMXSTP",	"MAX#_STEP",	h_ICHAR_f_int,	300,	"millibohr"},
{/*54*/  "IVIB",	"VIB#RATION",	h_ICHAR_handle,	0,	""},
{/*55*/  "IRECAL",	"EVA#L_HESS",	h_ICHAR_f_int,	-1,	"# of cyc."},
{/*56*/  "INTPROG",	"INTEGRAL#S",	h_ICHAR_handle,	1,	""},
{/*57*/  "IDISFD",	"FD_STEP#SIZE",	h_ICHAR_f_int,	0,	"10-4 bohr"},
{/*58*/  "IGRDFD",	"POI#NTS",	h_ICHAR_handle,	0,	""},
{/*59*/  "ICNTYP",	"CONT#RACTION",	h_ICHAR_handle,	1,	""},
{/*60*/  "ISYM",	"SYM#METRY",	h_ICHAR_handle,	2,	""},
/*       *alias,	*oppar,		ichar,		ideflt,	*units */
{/*61*/  "IBASIS",	"BAS#IS",	h_ICHAR_string,	0,	""},
{/*62*/  "IDFGHI",	"SPHER#ICAL",	h_ICHAR_handle,	1,	""},
{/*63*/  "IRESET",	"RESET#_FLAGS",	h_ICHAR_handle,	0,	""},
{/*64*/  "IPTORB",	"PER#T_ORB",	h_ICHAR_handle,	2,	""},
{/*65*/  "IGNBS1",	"GENBAS_1#",	h_ICHAR_f_int,	0,	""},
{/*66*/  "IGNBS2",	"GENBAS_2#",	h_ICHAR_f_int,	0,	""},
{/*67*/  "IGNBS3",	"GENBAS_3#",	h_ICHAR_f_int,	0,	""},
{/*68*/  "ICOORD",	"COORD#INATES",	h_ICHAR_handle,	3,	""},
{/*69*/  "ICKSYM",	"CHECK_SYM#",	h_ICHAR_handle,	1,	""},
{/*70*/  "ISCFPR",	"SCF_PR#INT",	h_ICHAR_f_int,	0,	""},
{/*71*/  "IECP",	"ECP#",		h_ICHAR_handle,	0,	""},
{/*72*/  "IRSTRT",	"RESTART#",	h_ICHAR_handle,	1,	""},
{/*73*/  "ITRAIN",	"TRA#NS_INV",	h_ICHAR_handle,	0,	""},
{/*74*/  "ISTABL",	"HFS#TABILITY",	h_ICHAR_handle,	0,	""},
{/*75*/  "ROTVEC",	"ROT_E#VEC",	h_ICHAR_f_int,	0,	""},
{/*76*/  "IBRTOL",	"BRUCK#_CONV",	h_ICHAR_f_int,	4,	"(tol)"},
{/*77*/  "IQRHFG",	"QRHF_G#ENERA",	h_ICHAR_string,	0,	""},
{/*78*/  "IUNITS",	"UNI#TS",	h_ICHAR_handle,	0,	""},
{/*79*/  "IFDGRP",	"FD_U#SEGROUP",	h_ICHAR_handle,	0,	""},
{/*80*/  "IFDPRJ",	"FD_P#ROJECT",	h_ICHAR_handle,	0,	""},
/*       *alias,	*oppar,		ichar,		ideflt,	*units */
/* 81    "IFDCAL",	"FD_C#ALTYPE",	h_ICHAR_handle,	0,	""*/
{/*81*/  "","",0,0,""},
{/*82*/  "IFDIRR",	"FD_IRR#EPS",	h_ICHAR_string,	0,	""},
{/*83*/  "IVTRAN",	"VTR#AN",	h_ICHAR_handle,	0,	""},
{/*84*/  "IHF2Fl",	"HF2#_FILE",	h_ICHAR_handle,	1,	""},
{/*85*/  "ISUBGP",	"SUBGRO#UP",	h_ICHAR_handle,	0,	""},
{/*86*/  "ISBXYZ",	"SUBGRP#AXIS",	h_ICHAR_handle,	0,	""},
{/*87*/  "IEXCIT",	"EXCITE#",	h_ICHAR_handle,	0,	""},
{/*88*/  "IZTACN",	"ZETA_CON#V",	h_ICHAR_f_int,	12,	"(tol)"},
{/*89*/  "IEXSYM",	"ESTATE_SYM#",	h_ICHAR_string,	0,	""},
{/*90*/  "ITREAT",	"TREA#T_PERT",	h_ICHAR_handle,	0,	""},
{/*91*/  "IEXPRP",	"ESTATE_PROP#",	h_ICHAR_handle,	0,	""},
{/*92*/  "IOPTCY",	"OPT_MAX#CYC",	h_ICHAR_f_int,	50,	""},
{/*93*/  "IABCDT",	"ABCDTYP#E",	h_ICHAR_handle,	0,	""},
{/*94*/  "IQRHFS",	"QRHF_S#PIN",	h_ICHAR_string,	0,	""},
{/*95*/  "IAOLAD",	"AO_LAD#DERS",	h_ICHAR_handle,	1,	""},
{/*96*/  "IFOCK",	"FOCK#",	h_ICHAR_handle,	1,	""},
{/*97*/  "IEXMXC",	"ESTATE_MAX#C",	h_ICHAR_f_int,	20,	""},
{/*98*/  "IEXTOL",	"ESTATE_TOL#",	h_ICHAR_f_int,	-1,	"(tol)"},
{/*99*/  "ITURBO",	"TURBOMOL#E",	h_ICHAR_handle,	0,	""},
{/*100*/ "IGABCD",	"GAMMA#_ABCD",	h_ICHAR_handle,	0,	""},
/*       *alias,	*oppar,		ichar,		ideflt,	*units */
{/*101*/ "IZTYPE",	"ZETA_TYP#E",	h_ICHAR_handle,	1,	""},
{/*102*/ "IZMAXC",	"ZETA_MAX#CYC",	h_ICHAR_f_int,	50,	""},
{/*103*/ "IRESRM",	"RESRAMAN#",	h_ICHAR_handle,	0,	""},
{/*104*/ "IPSI",	"PSI#",		h_ICHAR_handle,	0,	""},
{/*105*/ "IGEOPT",	"GEOM_OPT#",	h_ICHAR_handle,	0,	""},
{/*106*/ "IEXTRN",	"EXTERNAL#",	h_ICHAR_handle,	0,	""},
{/*107*/ "IHESUP",	"HESS_UPD#ATE",	h_ICHAR_handle,	0,	""},
{/*108*/ "IINHES",	"INIT_HESS#IAN",h_ICHAR_handle,	0,	""},
{/*109*/ "IEXTRP",	"EXTRAP#OLATE",	h_ICHAR_handle,	0,	""},
{/*110*/ "","",0,0,""},
{/*111*/ "","",0,0,""},
{/*112*/ "","",0,0,""},
{/*113*/ "","",0,0,""},
{/*114*/ "","",0,0,""},
{/*115*/ "","",0,0,""},
{/*116*/ "","",0,0,""},
{/*117*/ "","",0,0,""},
{/*118*/ "","",0,0,""},
{/*119*/ "","",0,0,""},
{/*120*/ "","",0,0,""},
/*       *alias,	*oppar,		ichar,		ideflt,	*units */
{/*121*/ "","",0,0,""},
{/*122*/ "","",0,0,""},
{/*123*/ "","",0,0,""},
{/*124*/ "","",0,0,""},
{/*125*/ "","",0,0,""},
{/*126*/ "","",0,0,""},
{/*127*/ "","",0,0,""},
{/*128*/ "","",0,0,""},
{/*129*/ "","",0,0,""},
{/*130*/ "","",0,0,""},
{/*131*/ "","",0,0,""},
{/*132*/ "","",0,0,""},
{/*133*/ "","",0,0,""},
{/*134*/ "","",0,0,""},
{/*135*/ "","",0,0,""},
{/*136*/ "","",0,0,""},
{/*137*/ "","",0,0,""},
{/*138*/ "","",0,0,""},
{/*139*/ "","",0,0,""},
{/*140*/ "","",0,0,""},
/*       *alias,	*oppar,		ichar,		ideflt,	*units */
{/*141*/ "","",0,0,""},
{/*142*/ "","",0,0,""},
{/*143*/ "","",0,0,""},
{/*144*/ "","",0,0,""},
{/*145*/ "","",0,0,""},
{/*146*/ "","",0,0,""},
{/*147*/ "","",0,0,""},
{/*148*/ "","",0,0,""},
{/*149*/ "","",0,0,""},
{/*150*/ "","",0,0,""},
{/*151*/ "","",0,0,""},
{/*152*/ "","",0,0,""},
{/*153*/ "","",0,0,""},
{/*154*/ "","",0,0,""},
{/*155*/ "","",0,0,""},
{/*156*/ "","",0,0,""},
{/*157*/ "","",0,0,""},
{/*158*/ "","",0,0,""},
{/*159*/ "","",0,0,""},
{/*160*/ "","",0,0,""},
/*       *alias,	*oppar,		ichar,		ideflt,	*units */
{/*161*/ "","",0,0,""},
{/*162*/ "","",0,0,""},
{/*163*/ "","",0,0,""},
{/*164*/ "","",0,0,""},
{/*165*/ "","",0,0,""},
{/*166*/ "","",0,0,""},
{/*167*/ "","",0,0,""},
{/*168*/ "","",0,0,""},
{/*169*/ "","",0,0,""},
{/*170*/ "","",0,0,""},
{/*171*/ "","",0,0,""},
{/*172*/ "","",0,0,""},
{/*173*/ "","",0,0,""},
{/*174*/ "","",0,0,""},
{/*175*/ "","",0,0,""},
{/*176*/ "","",0,0,""},
{/*177*/ "","",0,0,""},
{/*178*/ "","",0,0,""},
{/*179*/ "","",0,0,""},
{/*180*/ "","",0,0,""},
/*       *alias,	*oppar,		ichar,		ideflt,	*units */
{/*181*/ "","",0,0,""},
{/*182*/ "","",0,0,""},
{/*183*/ "","",0,0,""},
{/*184*/ "","",0,0,""},
{/*185*/ "","",0,0,""},
{/*186*/ "","",0,0,""},
{/*187*/ "","",0,0,""},
{/*188*/ "","",0,0,""},
{/*189*/ "","",0,0,""},
{/*190*/ "","",0,0,""},
{/*191*/ "","",0,0,""},
{/*192*/ "","",0,0,""},
{/*193*/ "","",0,0,""},
{/*194*/ "","",0,0,""},
{/*195*/ "","",0,0,""},
{/*196*/ "","",0,0,""},
{/*197*/ "","",0,0,""},
{/*198*/ "","",0,0,""},
{/*199*/ "","",0,0,""},
{/*200*/ "","",0,0,""},
/*       *alias,	*oppar,		ichar,		ideflt,	*units */
{/*201*/ "IEACLC",	"EA_CALC#",	h_ICHAR_handle,	0,	""},
{/*202*/ "IEASYM",	"EA_SYM#",	h_ICHAR_string,	0,	""},
{/*203*/ "ITDHF",	"TDHF#",	h_ICHAR_handle,	0,	""},
{/*204*/ "IFNCTL",	"FUNCT#IONAL",	h_ICHAR_handle,	4,	""},
{/*205*/ "IEOMCY",	"EOM_MAXCY#C",	h_ICHAR_f_int,	50,	"cycles"},
{/*206*/ "IEOMPR",	"EOMPROP#",	h_ICHAR_handle,	0,	""},
{/*207*/ "IABCDF",	"ABCDFULL#",	h_ICHAR_handle,	0,	""},
{/*208*/ "IINTOL",	"INTGRL_TOL#",	h_ICHAR_f_int,	14,	"(tol)"},
{/*209*/ "IDMPTY",	"DAMP_TYP#",	h_ICHAR_handle,	0,	""},
{/*210*/ "IDMPTL",	"DAMP_TOL#",	h_ICHAR_f_int,	10,	"x 0.01"},
{/*211*/ "ILSHA1",	"LSHF_A1#",	h_ICHAR_f_int,	0,	"x 0.01"},
{/*212*/ "ILSHB1",	"LSHF_B1#",	h_ICHAR_f_int,	0,	"x 0.01"},
{/*213*/ "IPOLYR",	"POLYRATE#",	h_ICHAR_handle,	0,	""},
{/*214*/ "IIPCLC",	"IP_CALC#",	h_ICHAR_handle,	0,	""},
{/*215*/ "IIPSYM",	"IP_SYM#",	h_ICHAR_string,	0,	""},
{/*216*/ "IPTYPE",	"IP_SEARCH#",	h_ICHAR_handle,	0,	""},
{/*217*/ "IEOM",	"EOMREF#",	h_ICHAR_handle,	0,	""},
{/*218*/ "ISOLV",	"SOLVEN#T",	h_ICHAR_f_int,	0,	""},
{/*219*/ "EETYPE",	"EE_SEARCH#",	h_ICHAR_handle,	0,	""},
{/*220*/ "IEOMPR",	"EOM_PRJCT#",	h_ICHAR_handle,	0,	""},
/*       *alias,	*oppar,		ichar,		ideflt,	*units */
{/*221*/ "INWVRT",	"NEWVRT#",	h_ICHAR_handle,	0,	""},
{/*222*/ "IABCD",	"HBARABCD#",	h_ICHAR_handle,	0,	""},
{/*223*/ "IABCI",	"HBARABCI#",	h_ICHAR_handle,	0,	""},
{/*224*/ "INT3EE",	"NT3EOMEE#",	h_ICHAR_handle,	0,	""},
{/*225*/ "INOREO",	"NOREORI#",	h_ICHAR_handle,	2,	""},
{/*226*/ "IEESYM",	"EE_SYM#",	h_ICHAR_string,	0,	""},
{/*227*/ "IKSPOT",	"KS_POT#",	h_ICHAR_handle,	0,	""},
{/*228*/ "IDIPC",	"DIP_CAL#C",	h_ICHAR_handle,	0,	""},
{/*229*/ "IDIPSY",	"DIP_SY#M",	h_ICHAR_string,	0,	""},
{/*230*/ "IDEAC",	"DEA_CAL#C",	h_ICHAR_handle,	0,	""},
{/*231*/ "IDEASY",	"DEA_SY#M",	h_ICHAR_string,	0,	""},
{/*232*/ "IPROG",	"PROGRAM#",	h_ICHAR_handle,	0,	""},
{/*233*/ "ICCR12",	"CCR12#",	h_ICHAR_handle,	0,	""},
{/*234*/ "IXEOMF",	"EOMXFIELD#",	h_ICHAR_f_int,	0,	"x 10-6"},
{/*235*/ "IYEOMF",	"EOMYFIELD#",	h_ICHAR_f_int,	0,	"x 10-6"},
{/*236*/ "IZEOMF",	"EOMZFIELD#",	h_ICHAR_f_int,	0,	"x 10-6"},
{/*237*/ "IINSF",	"INSERTF#",	h_ICHAR_handle,	0,	""},
{/*238*/ "IGRDCL",	"GRAD_CALC#",	h_ICHAR_handle,	0,	""},
{/*239*/ "IIMEM",	"IMEM#_SIZE",	h_ICHAR_f_int,	3000000,"Words"},
{/*240*/ "IMKRHF",	"MAKERHF#",	h_ICHAR_handle,	0,	""},
/*       *alias,	*oppar,		ichar,		ideflt,	*units */
{/*241*/ "IGLBMM",	"GLOBAL_MEM#",	h_ICHAR_f_int,	0,	"Words"},
{/*242*/ "IPRPNT",	"PRP_INT#S",	h_ICHAR_handle,	0,	""},
{/*243*/ "IFNOKP",	"FNO_KEEP#",	h_ICHAR_f_int,	0,	"percent"},
{/*244*/ "IFNOPT",	"FNO_POST#",	h_ICHAR_handle,	0,	""},
{/*245*/ "IFNOAC",	"FNO_ACT#IVE",	h_ICHAR_f_int,	0,	"percent"},
{/*246*/ "INAT",	"NATURAL#",	h_ICHAR_handle,	0,	""},
{/*247*/ "ICCSYM",	"ACC_SYM#",	h_ICHAR_string,	0,	""},
{/*248*/ "IUNO_R",	"UNO_REF#",	h_ICHAR_handle,	0,	""},
{/*249*/ "IUNO_C",	"UNO_CHARG#E",	h_ICHAR_f_int,	0,	""},
{/*250*/ "IUNO_M",	"UNO_MULT#",	h_ICHAR_f_int,	1,	""},
{/*251*/ "IRAMAN",	"RAMAN#",	h_ICHAR_handle,	0,	""},
{/*252*/ "IKUCH",	"KUCHARSKI#",	h_ICHAR_handle,	0,	""},
{/*253*/ "ISCF",	"SCF_TYPE#",	h_ICHAR_handle,	0,	""},
{/*254*/ "IDIRCT",	"DIRECT#",	h_ICHAR_handle,	0,	""},
{/*255*/ "BSNGST",	"SINGLE_STOR#E",h_ICHAR_handle,	0,	""},
{/*256*/ "","",0,0,""},
{/*257*/ "","",0,0,""},
{/*258*/ "","",0,0,""},
{/*259*/ "","",0,0,""},
{/*260*/ "","",0,0,""}
}; /* end ASV_nl[] definition */
#define MAX_ASVs 260

/******************************************************************************/

/* a tool for counting the NUMBER of strings in a NULL-string terminated
   char array (NOT including the NULL-string string) */
f_int count_strings(const char * handles[])
{
    f_int num = 0;
    while (*handles[num]) num++;
    return num;
}

/* a tool for updating an ASV of type handle based on a value token
   conversion */
void asv_update_handle(const f_int * index,
                       const char * value,
                       const char * handles[]
                      )
{
    f_int pos = 0;
    while ((strcasecmp(value,handles[pos])) && (*handles[pos])) pos++;
    if (!*handles[pos])
    {
     /* the value token is a number or junk */
        pos = atol(value);
        if (pos==0)
        {
            if (strcmp(value,"0"))
            {
                printf("\n     ERROR: The value token is not recognized.\n");
                F77_NAME(errex, ERREX)();
            }
        }
        else if ((pos<0) || (count_strings(handles)<=pos))
        {
            printf("\n     ERROR: %li is outside the acceptable interval "
                   "[0,%li]\n",pos,count_strings(handles)-1);
            F77_NAME(errex, ERREX)();
        }
    }
    f_flags.ioppar[*index] = pos;
    if (bVerbose_Update) printf("'%s' converts to '%s' -> %li\n",
                                value,
                                handles[f_flags.ioppar[*index]],
                                f_flags.ioppar[*index]);
    return;
}

/* this routine processes the ASVs of type handle and passes
   their handle arrays to asv_update_handle for updating */
void asv_handle_proc(const f_int * index, const char * value)
{
    switch (*index+1) /* Fortran-style indexing (to match h_IOPPAR_* )*/
    {

     /* You do not HAVE to have a separate case block for your ASV
        unless the list of handles is unique or requires special
        processing. If you find a handles[] array that matches your
        ASV, then PLEASE use it. Likewise, if you have to edit an
        existing array that belongs to more than one ASV, BE CAREFUL
        and create a new case block if necessary. */

     /* REMEMBER TO TERMINATE THE LISTS WITH NULL STRINGS! */

        case h_IOPPAR_calclevel:
        case h_IOPPAR_fno_keep:
        {
            const char *handles[] =
            {
      /* 0- 3*/ "SCF",		"MBPT(2)",	"MBPT(3)",	"SDQ-MBPT(4)",
      /* 4- 7*/ "MBPT(4)",	"LCCD",		"LCCSD",	"UCCSD(4)",
      /* 8-11*/ "CCD",		"UCC(4)",	"CCSD",		"CCSD[T]",
      /*12-15*/ "CCSD+TQ*",	"CCSDT-1",	"CCSDT-1b",	"CCSDT-2",
      /*16-19*/ "CCSDT-3",	"CCSDT-4",	"CCSDT",	"LCCSDT",
      /*20-23*/ "CCD+ST(CCD)",	"QCISD(T)",	"CCSD(T)",	"QCISD",
      /*24-27*/ "CID",		"CISD",		"QCISD(TQ)",	"CCSD(TQ)",
      /*28-31*/ "CCSD+TQ",	"CCSDT+Q*",	"CCSDT+Q",	"CC5SD(T)",
      /*32-35*/ "CCSD-T",	"CC3",		"CCSDT-T1T2",	"CCSDTQ-1",
      /*36-39*/ "CCSDTQF-1",	"CCSDTQ-2",	"CCSDTQ-3",	"CCSDTQ",
      /*40-43*/ "ACCSD",	"HFDFT",	"ACCSD(T)",	"CCSD(TQf)",
      /*44-47*/ "CCSDT(Qf)",	""
            };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_reference:
        {
            const char *handles[] = { "RHF",  "UHF",
                                      "ROHF", "TWODET", "ROHF-OS", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_props:
        {
            const char *handles[] =
            {
                "OFF",		/* 0 */
                "FIRST_ORDER",	/* 1 */
                "SECOND_ORDER",	/* 2 */
                "NMR",		/* 3 */
                "NMR_SWITCH",	/* 4 */
                "CONVENTIONAL",	/* 5 */
                "GIAO",		/* 6 */
                "TDHF",		/* 7 */
                "J_SO",		/* 8 */
                "J_FC",		/* 9 */
                "J_SD",		/* 10 */
                "EOM_NLO",	/* 11 */
                "GEERTSEN",	/* 12 */
                "JSC_ALL",	/* 13 */
                ""
            };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_density:
        {
            const char *handles[] = { "RELAXED", "RESPONSE", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_cc_extrapol:
        {
            const char *handles[] = { "STANDARD", "DIIS", "NOJACOBI",
                                      "OFF", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_scf_extrap:
        {
            const char *handles[] = { "NONE", "RPP", "QC",""};
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_incore:
        {
            const char *handles[] =
            {
                "NONE",		/* 0 */
                "NOABCD",	/* 1 */
                "T",		/* 2 */
                "=",		/* 3 (unused) */
                "NOABCI",	/* 4 */
                "NOABIJ",	/* 5 */
                "ALL",		/* 6 */
                ""
            };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_orbitals:
        {
            const char *handles[] = { "STANDARD", "SEMICANONICAL", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_guess:
        {
            const char *handles[] =
            {
                "MOREAD",	/* 0 */
                "CORE",		/* 1 */
                "NDDO",		/* 2 */
                "WALT_PRJDEN",	/* 3 */
                "READ_SO_MOS",	/* 4 */
                "READ_AO_MOS",	/* 5 */
                "MIN_BASIS",	/* 6 */
                "HUCKEL",	/* 7 */
                "OVERRIDE"      /* 8 */
                ""
            };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_geom_opt:
        {
            const char *handles[] =
            {
                "NONE",		/* 0 */
                "PARTIAL",   	/* 1 */
                "CART",		/* 2 */
                "FULL",		/* 3 */
                "RIC" ,         /* 4 */
                ""
            };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_external:
        {
            const char *handles[] =
            {
                "NONE",		/* 0 */
                "HYPERCHEM",	/* 1 */
                "MOLDEN",	/* 2 */
                ""
            };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_stp_siz_ctl:
        {
            const char *handles[] =
            {
                "TRUST_RADIUS",	/* 0 */
                "NORM",		/* 1 */
                "MAXIMUM",	/* 2 */
                ""
            };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_hess_update:
        {
            const char *handles[] =
            {
                "NONE",		/* 0 */
                "POWELL",	/* 1 */
                "BFGS",		/* 2 */
                "MS",		/* 3 */
                "BOFILL",	/* 4 */
                "PSB",		/* 5 */
                ""
            };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_init_hessian:
        {
            const char *handles[] =
            {
                "SPECIAL",	/* 0 */
                "FCMINT",	/* 1 */
                "MOPAC",	/* 2 */
                "EXACT",	/* 3 */
                ""
            };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_extrap:
        {
            const char *handles[] =
            {
                "NONE",		/* 0 */
                "GRADIENT",	/* 1 */
                "ENERGY",	/* 2 */
                "COMBO",	/* 3 */
                ""
            };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_opt_method:
        {
            const char *handles[] =
            {
                "AUTO",		/* 0 */
                "NR",		/* 1 */
                "RFA",		/* 2 */
                "MANR",		/* 3 */
                "EVFTS",	/* 4 */
                "MAEVFTS",	/* 5 */
                "IGTS",		/* 6 - used to be ENERONLY? */
                "QSD",		/* 7 */
                ""
            };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_negeval:
        {
            const char *handles[] = { "ABORT", "ABSVAL", "RFA", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_curvilinear:
        {
            const char *handles[] = { "OFF", "ON", "NO", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_vibration:
        {
            const char *handles[] = { "NO", "EXACT", "FINDIF_OLD", "FINDIF",
                                      "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_integrals:
        {
            const char *handles[] = { "ARGOS", "VMOL", "HERMIT", "CADPAC",
                                      "SEWARD", "GAMESS", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_points:
        {
            const char *handles[] = { "DOUBLE", "SINGLE", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_contraction:
        {
            const char *handles[] = { "SEGMENTED", "GENERAL", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_symmetry:
        {
            const char *handles[] = { "NONE", "OFF", "ON", "FULL",
                                      ""};
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_pert_orb:
        {
            const char *handles[] = { "STANDARD", "CANONICAL", "UNKNOWN", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_coordinates:
        {
            const char *handles[] = { "INTERNAL", "CARTESIAN", "XYZ2INT",
                                      "AUTO", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_check_sym:
        {
            const char *handles[] = { "NORMAL", "OVERRIDE", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_trans_inv:
        {
            const char *handles[] = { "USE", "IGNORE", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_hfstability:
        {
            const char *handles[] = { "OFF", "ON", "FOLLOW", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_units:
        {
            const char *handles[] = { "ANGSTROM", "BOHR", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_fd_usegroup:
        {
            const char *handles[] = { "FULL", "COMP", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_fd_project:
        {
            const char *handles[] = { "ON", "OFF", "" }; /* WHY GOD!?! */
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_vtran:
        {
            const char *handles[] = { "FULL/PARTIAL", "FULL", "PARTIAL", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_hf2_file:
        {
            const char *handles[] = { "SKIP", "USE", "SAVE", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_subgroup:
        {
            const char *handles[] =
            {
                "DEFAULT",	/* 0 */
                "C1",		/* 1 */
                "C2",		/* 2 */
                "CS",		/* 3 */
                "CI",		/* 4 */
                "C2V",		/* 5 */
                "C2H",		/* 6 */
                "D2",		/* 7 */
                "D2H",		/* 8 */
                ""
            };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_subgrpaxis:
        {
            const char *handles[] = { "X", "Y", "Z", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_excite:
        {
            const char *handles[] =
            {
                "NONE",		/* 0 */
                "TDA",		/* 1 */
                "RPA",		/* 2 */
                "EOMEE",	/* 3 */
                "-=EOMIP=-",	/* 4 - unused? */
                "CIS",		/* 5 */
                "CIS(D)",	/* 6 */
                "P-EOMEE",	/* 7 */
                "EOM-BWPT2",	/* 8 */
                "STEOM",	/* 9 */
                ""
            };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_treat_pert:
        {
            const char *handles[] = { "SIMULTANEOUS", "SEQUENTIAL", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_estate_prop:
        {
            const char *handles[] = { "OFF",       "EXPECTATION",
                                      "UNRELAXED", "DERIVATIVE",  "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_abcdtype:
        {
            const char *handles[] = { "STANDARD", "MULTIPASS", "AOBASIS", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_ao_ladders:
        {
            const char *handles[] = { "MULTIPASS", "SINGLEPASS", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_fock:
        {
            const char *handles[] = { "PK", "AO", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_gamma_abcd:
        {
            const char *handles[] = { "DISK", "DIRECT", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_zeta_type:
        {
            const char *handles[] = { "POPLE", "DIIS", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_ea_calc:
        {
            const char *handles[] =
            {
                "NONE",		/* 0 */
                "MBPT(2)",	/* 1 */
                "SO_DYSON",	/* 2 */
                "OVGF",		/* 3 */
                "P_EOMEA",	/* 4 */
                "EA_EOMCC",	/* 5 */
                "COMBO",	/* 6 */
                "OS_CCSD",	/* 7 */
                ""
            };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_functional:
        {
            const char *handles[] =
            {
                "XALPHA",	/* 0 */
                "BECKE",	/* 1 */
                "LYP",		/* 2 */
                "XLYP",		/* 3 */
                "BLYP",		/* 4 */
                ""
            };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_eomprop:
        {
            const char *handles[] = { "CILIKE", "LINEAR", "QUADRATIC", "COMBO",
                                      "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_abcdfull:
        {
            const char *handles[] = { "UNKNOWN", "ON", "OFF", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_damp_typ:
        {
            const char *handles[] = { "NONE", "DAVIDSON", "OTHER", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_ip_calc:
        {
            const char *handles[] =
            {
                "NONE",		/* 0 */
                "MBPT(2)",	/* 1 */
                "SO_DYSON",	/* 2 */
                "OVGF",		/* 3 */
                "P_EOMIP",	/* 4 */
                "IP_EOMCC",	/* 5 */
                "COMBO",	/* 6 */
                "OS_CCSD",	/* 7 */
                ""
            };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_ip_search:
        {
            const char *handles[] = { "VALENCE",  "LOWEST", "COREIP", "SHAKEUP",
                                      "KOOPMANS", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_eomref:
        {
            const char *handles[] = { "NONE", "CCSD", "MBPT(2)", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_ee_search:
        {
            const char *handles[] = { "LOWEST", "CORE", "LUMO", "HOMO", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_eom_prjct:
        {
            const char *handles[] = { "NO", "SEARCH_ONLY", "PRJCT_ALL",
                                      "PRJCT_NOISE", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_hbarabcd:
        case h_IOPPAR_hbarabci:
        {
            const char *handles[] = { "UNKNOWN", "OFF", "ON", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_nt3eomee:
        {
            const char *handles[] =
            {
                "NONE",		/* 0 */
                "NCCSDT-1",	/* 1 */
                "F-NCCSDT-1",	/* 2 */
                "NCCSDT-1b",	/* 3 */
                "F-NCCSDT-1b",	/* 4 */
                "CCSD(TPR)",	/* 5 */
                "F-CCSD(TPR)",	/* 6 */
                "CCSDR(T)",	/* 7 */
                "F-CCSDR(T)",	/* 8 */
                "CCSDR(3)",	/* 9 */
                "F-CCSDR(3)",	/* 10 */
                "NCCSDT-3",	/* 11 */
                "F-NCCSDT-3",	/* 12 */
                ""
            };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_ks_pot:
        {
            const char *handles[] =
            {
                "HF",		/* 0 */
                "XALPHA",	/* 1 */
                "LYP",		/* 2 */
                "XALYP",	/* 3 */
                "BLYP",		/* 4 */
                "LDA",		/* 5 */
                "B3LYP",	/* 6 */
                ""
            };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_dip_calc:
        case h_IOPPAR_dea_calc:
        {
            const char *handles[] = { "NONE", "TDA", "EOMCC", "STEOM",
                                      "OS_CCSD", "SS_STEOM", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_program:
        {
            const char *handles[] = { "DEFAULT", "NOT_USED", "ACES2", "MN_A3",
                                      "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_insertf:
        {
            const char *handles[] =
            {
                "OFF",		/* 0 */
                "SCF",		/* 1 */
                "TDA",		/* 2 */
                "CC",		/* 3 */
                "EOM",		/* 4 */
                "STEOM",	/* 5 */
                ""
            };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_grad_calc:
        {
            const char *handles[] = { "NONE", "ANALYTICAL", "NUMERICAL",""};
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_prp_ints:
        {
            const char *handles[] = { "PARTIAL", "FULL", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_scf_type:
        {
            const char *handles[] = { "HF", "KS", "HFDFT", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        case h_IOPPAR_noreori:
        {
            const char *handles[] = { "OFF", "ON", "AUTO", "" };
            asv_update_handle(index,value,handles);
            break;
        }

        default:
        {
         /* define default switching types */
            const char *one[] = { "ON",  "TRUE",   "YES",    "POSITIVE", "1",
                                  "ONE", "FIRST",  "SINGLE", "PRIMARY", "" };
            const char *two[] = { "TWO", "SECOND", "DOUBLE", "SECONDARY", "" };
            int pos;
            f_flags.ioppar[*index] = atol(value); /* initialize to an integer */
            pos = 0;
            while ((strcasecmp(value,one[pos])) && (*one[pos])) pos++;
            if (*one[pos])
            {
                f_flags.ioppar[*index] = 1;
                if (bVerbose_Update) printf("'%s' converts to 'ON' -> 1\n",
                                            value);
            }
            pos = 0; /* it is not efficient to continue scanning if we have
                        a match already, but hey, there is only one more... */
            while ((strcasecmp(value,two[pos])) && (*two[pos])) pos++;
            if (*two[pos])
            {
                f_flags.ioppar[*index] = 2;
                if (bVerbose_Update) printf("'%s' converts to 'TWO' -> 2\n",
                                            value);
            }
            if (!f_flags.ioppar[*index])
                if (bVerbose_Update) printf("'%s' converts to 'OFF' -> 0\n",
                                            value);
            break;
        }

    } /* end switch (*index+1) */
    return;
}

/******************************************************************************/

/* a tool for copying a full key definition into a string without the
   stub delimiter */
char *get_keyname(char *keyname, const char *keydef, const char *delim)
{
    size_t stublen = strcspn(keydef,delim);
    if (stublen)
    {
        strncpy(keyname,keydef,stublen);
        strcpy(keyname+stublen,keydef+stublen+1);
    }
    else
        strcpy(keyname,keydef);
    return keyname;
}

/******************************************************************************/

/* a tool for printing an ASV structure */
void
F77_NAME(asv_print,ASV_PRINT)
(const f_int * plASV)
{
    f_int index = *plASV-1;
    if ((0<*plASV) && (index<MAX_ASVs))
    {
        if (*ASV_nl[index].oppar)
        {
            char keyname[MAX_KEY_LENGTH];
            get_keyname(&keyname[0],ASV_nl[index].oppar,STUB_DELIM);
            printf("ASV[%3li] : ALIAS  = %s\n"
                   "           OPPAR  = %s (%s)\n"
                   "           ICHAR  = %li (type %s)\n"
                   "           IDEFLT = %li %s\n"
                   "       currently -> %li\n\n",
                   *plASV,
                   ASV_nl[index].alias,
                   ASV_nl[index].oppar, &keyname,
                   ASV_nl[index].ichar, ichar_types[ASV_nl[index].ichar],
                   ASV_nl[index].ideflt, ASV_nl[index].units,
                   f_flags.ioppar[index]
                  );
        }
    }
    return;
}

/******************************************************************************/

/* a tool for dumping the ASVs */
void F77_NAME(asv_dump,ASV_DUMP)()
{
    int iArg, i=0;
    printf(
"\n"
"          ASV#   ASV KEY DEFINITION =    CURRENT [   DEFAULT] UNITS\n"
"          ------------------------------------------------------------\n"
          );
    for (iArg=0;iArg<MAX_ASVs;iArg++)
    {
         if ((*ASV_nl[iArg].oppar) &&
             ((ASV_nl[iArg].ichar==h_ICHAR_handle) ||
              (ASV_nl[iArg].ichar==h_ICHAR_f_int)
             )
            )
         {
             char keyname[MAX_KEY_LENGTH];
             get_keyname(&keyname[0],ASV_nl[iArg].oppar,STUB_DELIM);
             printf("          %3i: %20s = %10li [%10li] %s\n",
                    iArg+1, &keyname, f_flags.ioppar[iArg],
                    ASV_nl[iArg].ideflt,ASV_nl[iArg].units
                   );
          /* print a blank line after every fourth printf */
             if (!(((i++)+1) & 3)) printf("\n");
         }
    }
    printf(
"          ------------------------------------------------------------\n"
"\n"
          );
    fflush(stdout);
    return;
}

/******************************************************************************/

/*
 * the key-value parser that converts a key token into an ASV index and:
 *  - updates the f_flags global structure with the value token (if type f_int),
 *  - passes the index and value token to a handle processor (if type handle),
 *  - copies the value token into the data location (if type string),
 *  - translates the value token into the data location
 *    (if type f_int array, double array, or double)
 */

/*
 * kvpair : a pointer to the whole "key=value" NULL TERMINATED string
 *
 *          IMPORTANT ! IMPORTANT ! IMPORTANT ! IMPORTANT ! IMPORTANT
 *          If the function is called by a Fortran routine, then it
 *          must append ACHAR(0) after the value token.
 *
 * index  : a pointer to the (FORTRAN!) index of the matched ASV
 *        = 0 if key does not match an ASV
 *        > 0 if key does     match an ASV
 *
 * data   : a pointer to the data output if the value type is string or array
 * size   : a pointer to the number of BYTES available at data
 *          (input)  the maximum size at data
 *          (output) the amount of data filled in to data
 */

void
F77_NAME(asv_update_kv,ASV_UPDATE_KV)
(const char * kvpair0, f_int * index, void * data, f_int * size)
{
    char *key, *value, kvpair[81];
    int NotKey;

 /* copy the input key-value token string to the kvpair buffer */
    strncpy(kvpair,kvpair0,81);
    if (kvpair[80])
    {
        printf("\n     WARNING: The key-value token string has been truncated "
                               "before processing.\n");
        kvpair[80]='\0';
    }

 /* point to the value token or be NULL */
    if (value=strpbrk(kvpair,"="))
    {
        *value='\0';
        value++;
        while (isspace(*value)) value++;
        if (*value)
        {
            char *s1 = value+strlen(value)-1;
            while (isspace(*s1)) s1--;
            *(s1+1)='\0';
        }
        else
            value='\0';
    }

 /* point to the key token or be NULL */
 /* If the first char is a '!', assume the bool to be false and go on. */
    key=kvpair;
    while (isspace(*key)) key++;
    NotKey = (*key=='!');
    if (NotKey)
    {
        if (value)
        {
         /* There had better not be a value token with a negation. */
            printf("\n     ERROR: Negation flag found with value token.\n"
                     "     key = '%s'; value = '%s'\n",key,value);
            F77_NAME(errex, ERREX)();
        }
        key++; while (isspace(*key)) key++;
    }
    if (*key)
    {
        char *s1 = key;
        while (*s1 && (!isspace(*s1))) s1++;
        *s1='\0';
    }
    else
        key='\0';

 /* convert the value token into an environment variable */
    if (value)
    {
        if (!strncmp(value,"${",2))
        {
            char * cz;
            const char *szNameSet = "abcdefghijklmnopqrstuvwxyz"
                                    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                    "1234567890_[]";
            if (cz=strpbrk(value,"}"))
                *cz='\0';
            else
            {
                printf("\n     ERROR: missing \"}\" in value token\n"
                         "            key = '%s'; value = '%s'\n",key,value);
                F77_NAME(errex, ERREX)();
            }
            value += 2; /* skip the "${" characters */
            while (isspace(*value)) value++;
            if (!*value)
                value='\0';
            else
            {
                *(value+strspn(value,szNameSet))='\0';
                if (cz=getenv(value))
                {
                    while (isspace(*cz)) cz++;
                    if (strlen(cz) > 80)
                    {
                        printf("\n     WARNING: The value for the environment "
                               "variable '%s'\n"
                                 "              has been truncated before "
                               "processing.\n",value);
                        *(cz+80)='\0';
                    }
                    if (*cz)
                    {
                        value=cz;
                        cz = value+strlen(value)-1;
                        while (isspace(*cz)) cz--;
                        *(cz+1)='\0';
                    }
                    else
                        value='\0';
                }
                else
                    value='\0';
            }
        }
    }

 /* initialize the ASV index to ERROR/UNKNOWN */
    *index = 0;

    if (key)
    {
     /* convert key token into an ASV index */
        f_int var = -1 + atol(key);

        if (var == -1)
        {
         /* match the key token to a known key alias or key stub */
         /* NOTE: These loops will scan all the ASVs and print an error
                  message if the token matches more than one ASV. */
            f_int lASV;
#ifdef PRINT_ALL
            printf("\n'%s' matches: ",key);
#endif /* PRINT_ALL */

         /* first, match any alias */
            for (lASV=0;lASV<MAX_ASVs;lASV++)
            {
                if (!strcasecmp(key,ASV_nl[lASV].alias))
                {
                    if (var==-1)
                    {
#ifdef PRINT_ALL
                        printf("%s\n",ASV_nl[lASV].alias);
#endif /* PRINT_ALL */
                        var=lASV;
                    }
                    else
                    {
                        printf("     ERROR: key token matched '%s' after "
                               "matching '%s'\n",
                               ASV_nl[lASV].alias,ASV_nl[var].alias
                              );
                        F77_NAME(errex, ERREX)();
                    }
                } /* end if key token matches ASV alias */
            } /* end loop over ASV aliases */

         /* now, match the stub */
            if (var==-1)
            {
                for (lASV=0;lASV<MAX_ASVs;lASV++)
                {
                    size_t stublen = strcspn(ASV_nl[lASV].oppar,STUB_DELIM);
                    if (stublen)
                    {
                        if (!strncasecmp(key,ASV_nl[lASV].oppar,stublen))
                        {
                            if (var==-1)
                            {
#ifdef PRINT_ALL
                                printf("%s\n",ASV_nl[lASV].oppar);
#endif /* PRINT_ALL */
                                var=lASV;
                            }
                            else
                            {
                                printf("     ERROR: key token matched '%s' "
                                       "after matching '%s'\n",
                                       ASV_nl[lASV].oppar,ASV_nl[var].oppar
                                      );
                                F77_NAME(errex, ERREX)();
                            }
                        } /* end if key token matches ASV stub */
                    } /* end if key definition contains a stub */
                } /* end loop over ASV stubs */
            } /* end if key token still unmatched */

        } /* end if the key token was not an integer (may still be unknown) */

        if ((var<0)||(MAX_ASVs<=var))
        {
         /* complain if the key token is unmatched */
            printf("     ERROR: '%s' does not match a known ASV!\n",key);
            F77_NAME(errex, ERREX)();
        }
        else
        {
         /* return the Fortran index */
            *index = var+1;

#ifdef PRINT_ALL
         /* print the structure for the requested ASV */
            F77_NAME(asv_print,ASV_PRINT)(index);
#endif /* PRINT_ALL */

         /* process the value token if available */
            if (value)
            {
                char keyname[MAX_KEY_LENGTH];
                get_keyname(&keyname[0],ASV_nl[var].oppar,STUB_DELIM);
                if (bVerbose_Update) printf("     Updating (%s,%s): ",
                                            keyname,value);
                switch (ASV_nl[var].ichar)
                {
                    case h_ICHAR_handle:
                    {
                        asv_handle_proc(&var,value);
                        break;
                    }
                    case h_ICHAR_string:
                    {
                     /* measure the byte length of the value token w/o NULL */
                        size_t length = strlen(value);
                        size_t Bytes  = sizeof(*value)*length;
                        if (*size>=(f_int)Bytes)
                        {
                            memmove(data,value,Bytes);
                            *size=(f_int)Bytes;
                            if (bVerbose_Update)
                                printf("String copy succeeded.\n");
                        }
                        else
                        {
                            printf("String copy failed, only %li characters"
                                   " copied.\n",*size);
                            F77_NAME(errex, ERREX)();
                            memmove(data,value,(size_t)*size);
                        }
                        break;
                    }
                    case h_ICHAR_f_int:
                    {
                     /* f_int tmp = atol(value); */
                        f_int tmp;
                        long lTmp; sscanf(value,"%li",&lTmp); tmp = lTmp;
                        if ((!tmp) && strcmp(value,"0") && strcmp(value,"0x0"))
                        {
                            printf("\n     ERROR: The value token is not "
                                   "recognized.\n");
                            F77_NAME(errex, ERREX)();
                        }
                     /* This block is mainly for processing units. If any
                        other units are added besides data sizes, then
                        PLEASE move this section to a separate function. */
                        switch (var+1)
                        {
                            char * czUnit;
                            case h_IOPPAR_mem:
                            case h_IOPPAR_imem:
                            case h_IOPPAR_file_recsiz:
                            case h_IOPPAR_global_mem:
                             /* scale the value by any supplied unit */
                                czUnit = value;
                                while ((!isalpha(*czUnit)) &&
                                       (*czUnit!='\0')        ) czUnit++;
                                if (*czUnit)
                                {
                                 /* allow for '1.5gb' */
                                    double dTmp = atof(value);
                                    switch (*czUnit)
                                    {
                                        case 'B': case 'b':
                                            dTmp /= sizeof(f_int);
                                            break;
                                        case 'W': case 'w':
                                            break;
                                        default:
                                         /* switch on a multiplier */
                                            switch (*czUnit)
                                            {
                                                case 'M': case 'm':
                                                    dTmp *= 1048576;    break;
                                                case 'G': case 'g':
                                                    dTmp *= 1073741824; break;
                                                case 'K': case 'k':
                                                    dTmp *= 1024;       break;
                                                default:
                                                    printf("ERROR: '%s' is not "
                                                           "a valid unit.\n",
                                                           czUnit);
                                                    F77_NAME(errex, ERREX)();
                                                    break;
                                            }
                                         /* switch on a unit */
                                            switch (*(czUnit+1))
                                            {
                                                case 'B': case 'b':
                                                    dTmp /= sizeof(f_int);
                                                    break;
                                                case 'W': case 'w':
                                                    break;
                                                default:
                                                    printf("ERROR: '%s' is not "
                                                           "a valid unit.\n",
                                                           czUnit);
                                                    F77_NAME(errex, ERREX)();
                                                    break;
                                            }
                                            break;
                                    } /* end if no multiplier */
                                    tmp = (f_int)dTmp;
                                } /* end if units found */
                                break;
                            default:
                                break;
                        } /* end ASV switch */
                        if (bVerbose_Update)
                            printf("%li -> ",f_flags.ioppar[var]);
                        f_flags.ioppar[var]=tmp;
                        if (bVerbose_Update)
                            printf("%li\n",  f_flags.ioppar[var]);
                        break;
                    }
                    case h_ICHAR_double:
                    {
                        double dTmp; char *s1;
                        size_t Bytes = sizeof(dTmp);
                        if (*size>=(f_int)Bytes)
                        {
                            if      (s1=strpbrk(value,"d")) *s1='e';
                            else if (s1=strpbrk(value,"D")) *s1='E';
                            dTmp = atof(value);
                            memmove(data,&dTmp,Bytes);
                            *size=(f_int)Bytes;
                            if (bVerbose_Update) printf("Returning %g\n",dTmp);
                        }
                        else
                        {
                            printf("ERROR: Type double expected of %i Bytes, "
                                   "but destination is only %li Bytes.\n",
                                   Bytes,*size
                                  );
                            F77_NAME(errex, ERREX)();
                            *size=0;
                        }
                        break;
                    }
                    case h_ICHAR_f_int_array:
                    {
                        printf("There is no array processing routine yet.\n");
                        break;
                    }
                    case h_ICHAR_double_array:
                    {
                        printf("There is no array processing routine yet.\n");
                        break;
                    }
                    default:
                    {
                        printf("ERROR: An unknown type has been entered!\n");
                        F77_NAME(errex, ERREX)();
                        break;
                    }
                } /* end ichar switch */
            }
            else /* no value token */
            {
                char keyname[MAX_KEY_LENGTH];
                get_keyname(&keyname[0],ASV_nl[var].oppar,STUB_DELIM);
                if (ASV_nl[var].ichar==h_ICHAR_handle)
                {
                    if (NotKey)
                    {
                        if (bVerbose_Update)
                            printf("     Updating (%s,OFF): ",keyname);
                        asv_handle_proc(&var,"OFF");
                    }
                    else
                    {
                        if (bVerbose_Update)
                            printf("     Updating (%s,ON): ",keyname);
                        asv_handle_proc(&var,"ON");
                    }
                }
                else
                {
                    printf("     ERROR: %s (%s) requires a value!\n",
                           key,keyname);
                    F77_NAME(errex, ERREX)();
                }
            } /* end if (value) */
        } /* end if ((var<0)||(MAX_ASVs<=var)) */

    }
    else /* key==NULL */
    {
     /* initialize the f_flags global structure */
        int iASV;
        for (iASV=0;iASV<MAX_ASVs;iASV++)
             f_flags.ioppar[iASV] = ASV_nl[iASV].ideflt;
        if (bVerbose_Update)
            printf("     The ACES State Variables were initialized.\n");
    }

    fflush(stdout);
    return;
}

/******************************************************************************/

void F77_NAME(init_flags, INIT_FLAGS) ()
{
   /* initialize the f_flags global structure */
   int iASV;

   for (iASV=0;iASV<MAX_ASVs;iASV++)
      f_flags.ioppar[iASV] = ASV_nl[iASV].ideflt;
}

int cli_prog(int argc, char *argv[])
/*int main(int argc, char *argv[])*/
{
    f_int iArg = 0;
    char szJunk[32]; f_int lSize = sizeof(szJunk);

 /* print a nice header */
    printf(
"\n"
"                      ACES STATE VARIABLE REGISTRATION LOG\n"
"     ----------------------------------------------------------------------\n"
          );

 /* initialize the f_flags global data */
    F77_NAME(asv_update_kv,ASV_UPDATE_KV)("",&iArg,&szJunk,&lSize);

 /* loop over arguments */
    for (iArg=1;iArg<argc;iArg++)
    {
         f_int var;
         F77_NAME(asv_update_kv,ASV_UPDATE_KV)(argv[iArg],&var,&szJunk,&lSize);
         /*
         if (var) printf("Update successful for ASV #%li\n\n",var);
         else     printf("Update failed for key-value '%s'\n\n",argv[iArg]);
         */
    }

 /* print a nice footer */
    printf(
"     ----------------------------------------------------------------------\n"
"\n"
          );

 /* dump the ASVs */
    F77_NAME(asv_dump,ASV_DUMP)();

    return 0;
}
