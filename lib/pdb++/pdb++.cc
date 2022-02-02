//
//	Copyright (c) 1992 The Regents of the University of California.
//	All rights reserved.
//
//	Redistribution and use in source and binary forms are permitted
//	provided that the above copyright notice and this paragraph are
//	duplicated in all such forms and that any documentation,
//	advertising materials, and other materials related to such
//	distribution and use acknowledge that the software was developed
//	by the University of California, San Francisco.  The name of the
//	University may not be used to endorse or promote products derived
//	from this software without specific prior written permission.
//	THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
//	IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
//	WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
//	$Id: pdb++.cc,v 1.6 1994/12/13 22:41:17 gregc Exp $
//

#include "pdb++.h"

extern "C" {
#include <string.h>
}

void
PDB::type(RecordType t)
{
	if (t == UNKNOWN) {
		// optimize default case (skip memset())
		rType = t;
		unknown.junk[0] = '\0';
		return;
	}
	memset(this, 0, sizeof *this);
	rType = t;
	switch (t) {
	default:
		break;
	case ATOM:
		atom.occupancy = 1.0;
		break;
	}
}

int
PDB::byteCmp(const PDB &l, const PDB &r)
{
	return memcmp(&l, &r, sizeof (PDB));
}
