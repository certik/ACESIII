#ifndef _F77_NAME_H_
#define _F77_NAME_H_

#ifdef C_UPPER
#      define F77_NAME(name,NAME) NAME
#else
#   ifdef C_SUFFIX
#      define F77_NAME(name,NAME) name##_
#   else
#      define F77_NAME(name,NAME) name
#   endif /* C_SUFFIX */
#endif /* C_UPPER */

#ifdef CB_UPPER
#      define F77_CB_NAME(name,NAME) NAME
#else
#   ifdef CB_SUFFIX
#      define F77_CB_NAME(name,NAME) name##_
#   else
#      define F77_CB_NAME(name,NAME) name
#   endif /* CB_SUFFIX */
#endif /* CB_UPPER */

#endif /* _F77_NAME_H_ */
