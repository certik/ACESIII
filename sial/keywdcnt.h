/******************************
keywdcnt.h
Lei Wang 
wang@qtp.ufl.edu
May 2004
******************************/

#define WILDCARD_INDICATOR 90909

#ifndef _keywdcnt_h
#define _keywdcnt_h

const int max_input_size=   5000000;
const int max_output_size=  1000000;
const int max_id_size=      1000000;
const int maxidlen=128;
const int operat_maxnumber =162;
const int noperat_keywords=108;

//define the max number of operators

enum operators {norb, nocc, nvirt, bocc,eocc,bvirt, evirt,
                naocc, nbocc, navirt, nbvirt, baocc, bbocc, 
                eaocc, ebocc, bavirt, bbvirt, eavirt, ebvirt,
                noccorb, nvirtorb, boccorb,eoccorb,bvirtorb, evirtorb,
                naoccorb, nboccorb, navirtorb, nbvirtorb, baoccorb, bboccorb,
                eaoccorb, eboccorb, bavirtorb, bbvirtorb, eavirtorb, ebvirtorb,
                cc_iter, cc_hist, cc_beg, scf_iter, scf_hist, scf_beg,
                natoms,
                itrips, itripe, ihess1, ihess2, jhess1, jhess2,
                subb, sube, sip_sub_segsize, sip_sub_occ_segsize,
                sip_sub_virt_segsize, sip_sub_ao_segsize,
                indexx, aoindex, moindex, moaindex, mobindex, laindex,
                subindex, window_index, served, temp, 
                distributed, staticc, local, userinstr, scalar,
                sial, endsial, proc, endproc,
                call, doo, enddo, iff, computeint,
                setindex, put, putadd, collective, execute, get,
                prepare, prepareadd, request, destroye, writee,
                create, deletee, allocatee, deallocatee, endiff, elsee,
                cycle, exitt,pardo, endpardo, also, returnsil,
                prequest, where, of_kwd, in_kwd,
                addassign, subassign, 
                multassign, divassign, 
				eq,ge, le, gt, lt, ne,
                assign, add, sub, mult, divide, tensor,
                lparen, rparen, lbracket, rbracket, comma, returnn,
                id, fconst, iconst, symbolic_const, jump, start,
                contract, sum, diff,
                eqf, gef, lef, gtf, ltf, nef,
				addf, subf, multf, dividef,
                idf, divide_array, error};

const char* const keywords[operat_maxnumber]=
                {"norb", "nocc", "nvirt", "bocc", "eocc", "bvirt", "evirt",
                "naocc", "nbocc", "navirt", "nbvirt", "baocc", "bbocc", 
                "eaocc", "ebocc", "bavirt", "bbvirt", "eavirt", "ebvirt",
                "noccorb", "nvirtorb", "boccorb","eoccorb",
                "bvirtorb", "evirtorb", "naoccorb", "nboccorb",
                "navirtorb", "nbvirtorb", "baoccorb", "bboccorb",
                "eaoccorb", "eboccorb", "bavirtorb", "bbvirtorb", 
                "eavirtorb", "ebvirtorb", 
                "cc_iter","cc_hist","cc_beg","scf_iter","scf_hist","scf_beg",
                "natoms",
                "itrips", "itripe", "ihess1", "ihess2", "jhess1", "jhess2",
                "subb", "sube", "sip_sub_segsize", "sip_sub_occ_segsize",
                "sip_sub_virt_segsize", "sip_sub_ao_segsize",
                "index", "aoindex", "moindex", "moaindex", "mobindex", "laindex",
                "subindex", "window_index", "served", "temp",
                "distributed", "static", "local", "userinstr", "scalar",
                "sial", "endsial", "proc", "endproc",
                "call", "do", "enddo", "if", "compute_integrals",
                "setindex", "put", "put", "collective", "execute", "get",
                "prepare", "prepare", "request", "destroy", "write",
                "create", "delete", "allocate", "deallocate", "endif", "else",
                "cycle", "exit", "pardo", "endpardo", "also", "return",
                "prequest", "where", "of", "in",
                "+=", "-=", "*=", "/=", 
				"==", ">=", "<=", ">", "<", "!=",
                "=", "+", "-", "*", "/", "^",
                "(", ")", "[", "]", ",", "\n",
                "id", "fconst", "iconst","symbolic_const", "jump", "start",
                "contract", "sum", "diff",
                "==", ">=", "<=", ">", "<", "!=",
				"+", "-", "*", "/",
                "idf", "/", "error"};
#endif
