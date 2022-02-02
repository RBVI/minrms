//
//	Copyright (c) 1989,1992 The Regents of the University of California.
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
//	$Id: pdb_chars.cc,v 1.8 95/02/28 14:09:52 gregc Exp $
//
//	subroutine for writing PDB format files
//

# include	"pdb++.h"
extern "C" {
# include	<ctype.h>
}
extern "C" int	sprintf(char *, const char *, ...);

static char const * const pdbRecordFormat[PDB::NUM_TYPES] = {
#include "write_format.i"
};

static char const * const pdbrun5[PDB::NUM_USER] = {
#include "pdbrun5_write.i"
};

static char const * const pdbrun6[PDB::NUM_USER] = {
#include "pdbrun6_write.i"
};

const char *
PDB::chars(void) const
{
	static char	buf[BufLen];
	const char	*fmt;
	const Sheet	*sh;
	const Residue	*shr0, *shr1, *sha0, *sha1;
	int		count;

	// convert C structure to pdb record

	if (rType < USER_PDBRUN)
		fmt = pdbRecordFormat[rType];
	else if (pdbrunOutputVersion < 6)
		fmt = pdbrun5[rType - USER_PDBRUN];
	else
		fmt = pdbrun6[rType - USER_PDBRUN];
	switch (rType) {

	case UNKNOWN:
		count = sprintf(buf, fmt, unknown.junk);
		break;

	case AGGRGT:
		count = sprintf(buf, fmt, aggrgt.serialNum,
			aggrgt.numComponents, aggrgt.cmpontSerialNums[0],
			aggrgt.cmpontSerialNums[1],
			aggrgt.cmpontSerialNums[2],
			aggrgt.cmpontSerialNums[3],
			aggrgt.cmpontSerialNums[4],
			aggrgt.cmpontSerialNums[5],
			aggrgt.cmpontSerialNums[6],
			aggrgt.cmpontSerialNums[7],
			aggrgt.cmpontSerialNums[8],
			aggrgt.cmpontSerialNums[9],
			aggrgt.cmpontSerialNums[10],
			aggrgt.cmpontSerialNums[11],
			aggrgt.cmpontSerialNums[12],
			aggrgt.cmpontSerialNums[13]);
		break;

	case AGRDES:
	case CMPDES:
	case FTNOTE:
	case MTXDES:
	case REMARK:
	case SYMDES:
		count = sprintf(buf, fmt, agrdes.num, agrdes.text);
		break;

	case ANISOU:
	case SIGUIJ:
		count = sprintf(buf, fmt, anisou.serialNum, anisou.name,
			anisou.altLoc, anisou.residue.name,
			anisou.residue.chainId, anisou.residue.seqNum,
			anisou.residue.insertCode, anisou.u[0], anisou.u[1],
			anisou.u[2], anisou.u[3], anisou.u[4], anisou.u[5]);
		break;

	case ATOM:
	case HETATM:
	case SIGATM:
		count = sprintf(buf, fmt, atom.serialNum, atom.name,
			atom.altLoc, atom.residue.name, atom.residue.chainId,
			atom.residue.seqNum, atom.residue.insertCode,
			atom.xyz[0], atom.xyz[1], atom.xyz[2], atom.occupancy,
			atom.tempFactor, atom.ftnoteNum);
		break;

	case AUTHOR:
	case COMPND:
	case JRNL:
	case SOURCE:
	case EXPDTA:
		count = sprintf(buf, fmt, author.continuation, author.data);
		break;

	case CONECT:
		count = sprintf(buf, fmt, conect.serialNum,
			conect.covalent[0], conect.covalent[1],
			conect.covalent[2], conect.covalent[3],
			conect.bonds[0].hydrogen[0],
			conect.bonds[0].hydrogen[1], conect.bonds[0].salt,
			conect.bonds[1].hydrogen[0],
			conect.bonds[1].hydrogen[1], conect.bonds[1].salt);
		break;

	case CMPONT:
		count = sprintf(buf, fmt, cmpont.seqNum,
			cmpont.residues[0].name, cmpont.residues[0].chainId,
			cmpont.residues[0].seqNum,
			cmpont.residues[0].insertCode,
			cmpont.residues[1].name, cmpont.residues[1].chainId,
			cmpont.residues[1].seqNum,
			cmpont.residues[1].insertCode);
		break;

	case CRYST1:
		count = sprintf(buf, fmt, cryst1.a, cryst1.b, cryst1.c,
			cryst1.alpha, cryst1.beta, cryst1.gamma,
			cryst1.spaceGrp, cryst1.z);
		break;

	case END:
	case ENDMDL:
		count = sprintf(buf, fmt);
		break;

	case FORMUL:
		count = sprintf(buf, fmt, formul.component, formul.hetId,
			formul.continuation, formul.exclude, formul.formula);
		break;

	case HEADER:
		count = sprintf(buf, fmt, header.classification,
			header.timestamp, header.type, header.id);
		break;

	case HELIX:
		count = sprintf(buf, fmt, helix.serialNum, helix.id,
			helix.residues[0].name, helix.residues[0].chainId,
			helix.residues[0].seqNum,
			helix.residues[0].insertCode, helix.residues[1].name,
			helix.residues[1].chainId, helix.residues[1].seqNum,
			helix.residues[1].insertCode, helix.type,
			helix.comment);
		break;

	case HET:
		count = sprintf(buf, fmt, het.hetGrp.name,
			het.hetGrp.chainId, het.hetGrp.seqNum,
			het.hetGrp.insertCode, het.numAtoms, het.text);
		break;

	case MASTER:
		count = sprintf(buf, fmt, master.numRemark, master.numFtnote,
			master.numHet, master.numHelix, master.numSheet,
			master.numTurn, master.numSite, master.numTransform,
			master.numCoordinate, master.numTer,
			master.numConect, master.numSeqres);
		break;

	case MODEL:
		count = sprintf(buf, fmt, model.num);
		break;

	case MTRIX:
		count = sprintf(buf, fmt, mtrix.rowNum, mtrix.serialNum,
			mtrix.m1, mtrix.m2, mtrix.m3, mtrix.v, mtrix.given);
		break;

	case OBSLTE:
		count = sprintf(buf, fmt, obslte.continuation, obslte.timestamp,
			obslte.oldId, obslte.idMap[0], obslte.idMap[1],
			obslte.idMap[2], obslte.idMap[3], obslte.idMap[4],
			obslte.idMap[2], obslte.idMap[6], obslte.idMap[7]);
		break;

	case ORIGX:
		count = sprintf(buf, fmt, origx.rowNum, origx.o1, origx.o2,
			origx.o3, origx.t);
		break;

	case REVDAT:
		count = sprintf(buf, fmt, revdat.modification,
			revdat.continuation, revdat.timestamp, revdat.id,
			revdat.modType, revdat.corrections);
		break;

	case SCALE:
		count = sprintf(buf, fmt, scale.rowNum, scale.s1, scale.s2,
			scale.s3, scale.u);
		break;

	case SEQRES:
		count = sprintf(buf, fmt, seqres.serialNum, seqres.chainId,
			seqres.count, seqres.names[0], seqres.names[1],
			seqres.names[2], seqres.names[3], seqres.names[4],
			seqres.names[5], seqres.names[6], seqres.names[7],
			seqres.names[8], seqres.names[9], seqres.names[10],
			seqres.names[11], seqres.names[12]);
		break;

	case SHEET:
		sh = &sheet;
		shr0 = &sh->residues[0];
		shr1 = &sh->residues[1];
		sha0 = &sh->atoms[0].residue;
		sha1 = &sh->atoms[1].residue;
		count = sprintf(buf, fmt, sh->strandNum, sh->id, sh->count,
			shr0->name, shr0->chainId, shr0->seqNum,
			shr0->insertCode, shr1->name, shr1->chainId,
			shr1->seqNum, shr1->insertCode, sh->sense,
			sh->atoms[0].name, sha0->name, sha0->chainId,
			sha0->seqNum, sha0->insertCode, sh->atoms[1].name,
			sha1->name, sha1->chainId, sha1->seqNum,
			sha1->insertCode);
		break;

	case SITE:
		shr0 = &site.residues[0];
		shr1 = &site.residues[1];
		sha0 = &site.residues[2];
		sha1 = &site.residues[3];
		count = sprintf(buf, fmt, site.seqNum, site.id, site.count,
			shr0->name, shr0->chainId, shr0->seqNum,
			shr0->insertCode,
			shr1->name, shr1->chainId, shr1->seqNum,
			shr1->insertCode,
			sha0->name, sha0->chainId, sha0->seqNum,
			sha0->insertCode,
			sha1->name, sha1->chainId, sha1->seqNum,
			sha1->insertCode);
		break;

	case SPRSDE:
		count = sprintf(buf, fmt, sprsde.continuation, sprsde.timestamp,
			sprsde.id, sprsde.supersede[0], sprsde.supersede[1],
			sprsde.supersede[2], sprsde.supersede[3],
			sprsde.supersede[4], sprsde.supersede[5],
			sprsde.supersede[6], sprsde.supersede[7]);
		break;

	case SSBOND:
		count = sprintf(buf, fmt, ssbond.seqNum,
			ssbond.residues[0].name, ssbond.residues[0].chainId,
			ssbond.residues[0].seqNum,
			ssbond.residues[0].insertCode,
			ssbond.residues[1].name, ssbond.residues[1].chainId,
			ssbond.residues[1].seqNum,
			ssbond.residues[1].insertCode, ssbond.comment);
		break;

	case SYMOP:
		count = sprintf(buf, fmt, symop.rowNum, symop.serialNum,
			symop.s1, symop.s2, symop.s3, symop.t);
		break;

	case TER:
		count = sprintf(buf, fmt, ter.serialNum, ter.residue.name,
			ter.residue.chainId, ter.residue.seqNum,
			ter.residue.insertCode);
		break;

	case TRNSFM:
		count = sprintf(buf, fmt, trnsfm.resultSerialNum,
			trnsfm.applySerialNum, trnsfm.sourceSerialNum);
		break;

	case TURN:
		count = sprintf(buf, fmt, turn.seqNum, turn.id,
			turn.residues[0].name, turn.residues[0].chainId,
			turn.residues[0].seqNum, turn.residues[0].insertCode,
			turn.residues[1].name, turn.residues[1].chainId,
			turn.residues[1].seqNum, turn.residues[1].insertCode,
			turn.comment);
		break;

	case TVECT:
		count = sprintf(buf, fmt, tvect.serialNum, tvect.t1, tvect.t2,
			tvect.t3, tvect.comment);
		break;

	case USER:
		count = sprintf(buf, fmt, user.subtype, user.text);
		break;

	case USER_PDBRUN:
		count = sprintf(buf, fmt, userPdbrun.version);
		pdbrunOutputVersion = userPdbrun.version;
		break;

	case USER_EYEPOS:
		count = sprintf(buf, fmt, userEyePos.xyz[0], userEyePos.xyz[1],
			userEyePos.xyz[2]);
		break;

	case USER_ATPOS:
		count = sprintf(buf, fmt, userAtPos.xyz[0], userAtPos.xyz[1],
			userAtPos.xyz[2]);
		break;

	case USER_WINDOW:
		count = sprintf(buf, fmt, userWindow.left, userWindow.right,
			userWindow.bottom, userWindow.top, userWindow.hither,
			userWindow.yon);
		break;

	case USER_FOCUS:
		count = sprintf(buf, fmt, userFocus.focus);
		break;

	case USER_VIEWPORT:
		count = sprintf(buf, fmt, userViewport.xmin, userViewport.xmax,
			userViewport.ymin, userViewport.ymax);
		break;

	case USER_BGCOLOR:
		if (pdbrunOutputVersion < 6)
			count = ::sprintf(buf, fmt, userBgColor.rgb[0],
				userBgColor.rgb[1], userBgColor.rgb[2]);
		else
			count = sprintf(buf, fmt, userBgColor.rgb[0],
				userBgColor.rgb[1], userBgColor.rgb[2]);
		break;

	case USER_ANGLE:
		if (pdbrunOutputVersion < 6)
			count = sprintf(buf, fmt, userAngle.which,
				userAngle.atom0, userAngle.atom1,
				userAngle.atom2, userAngle.atom3,
				userAngle.angle);
		else
			count = sprintf(buf, fmt, userAngle.atom0,
				userAngle.atom1, userAngle.atom2,
				userAngle.atom3, userAngle.angle);
		break;

	case USER_DISTANCE:
		if (pdbrunOutputVersion < 6)
			count = sprintf(buf, fmt, userDistance.which,
				userDistance.atom0, userDistance.atom1,
				userDistance.distance);
		else
			count = sprintf(buf, fmt, userDistance.atom0,
				userDistance.atom1, userDistance.distance);
		break;

	case USER_FILE:
		if (pdbrunOutputVersion < 6)
			count = sprintf(buf, fmt, userFile.filename);
		else
			count = sprintf(buf, fmt, userFile.model,
							userFile.filename);
		break;

	case USER_MARKNAME:
		if (pdbrunOutputVersion < 6)
			count = 0;
		else
			count = sprintf(buf, fmt, userMarkname.markname);
		break;

	case USER_MARK:
		if (pdbrunOutputVersion < 6)
			count = 0;
		else
			count = sprintf(buf, fmt, userMark.markname);
		break;

	case USER_CNAME:
		if (pdbrunOutputVersion < 6)
			count = ::sprintf(buf, fmt, userCName.name,
				userCName.rgb[0], userCName.rgb[1],
				userCName.rgb[2]);
		else
			count = sprintf(buf, fmt, userCName.rgb[0],
				userCName.rgb[1], userCName.rgb[2],
				userCName.name);
		break;

	case USER_COLOR:
		if (pdbrunOutputVersion < 6)
			count = ::sprintf(buf, fmt, userColor.spec,
				userColor.rgb[0], userColor.rgb[1],
				userColor.rgb[2]);
		else
			count = sprintf(buf, fmt, userColor.rgb[0],
				userColor.rgb[1], userColor.rgb[2],
				userColor.spec);
		break;

	case USER_RADIUS:
		count = sprintf(buf, fmt, userRadius.radius);
		break;

	case USER_OBJECT:
		count = sprintf(buf, fmt, userObject.model);
		break;

	case USER_ENDOBJ:
		count = sprintf(buf, fmt, userEndObj.model);
		break;

	case USER_CHAIN:
		if (pdbrunOutputVersion < 6)
			count = ::sprintf(buf, fmt, userChain.atom0,
				userChain.atom1);
		else
			count = sprintf(buf, fmt, userChain.atom0,
				userChain.atom1);
		break;

	case USER_GFX_BEGIN:
		if (pdbrunOutputVersion < 6)
			count = 0;
		else if (userGfxBegin.primitive == GFX_UNKNOWN)
			count = sprintf(buf, fmt, userGfxBegin.unknown);
		else
			count = sprintf(buf, fmt,
					gfxChars(userGfxBegin.primitive));
		break;

	case USER_GFX_END:
		if (pdbrunOutputVersion < 6)
			count = 0;
		else
			count = sprintf(buf, fmt);
		break;

	case USER_GFX_COLOR:
		if (pdbrunOutputVersion < 6)
			count = ::sprintf(buf, fmt, userGfxColor.spec,
				userGfxColor.rgb[0], userGfxColor.rgb[1],
				userGfxColor.rgb[2]);
		else
			count = sprintf(buf, fmt, userGfxColor.rgb[0],
				userGfxColor.rgb[1], userGfxColor.rgb[2],
				userGfxColor.spec);
		break;

	case USER_GFX_NORMAL:
		if (pdbrunOutputVersion < 6)
			count = 0;
		else
			count = sprintf(buf, fmt, userGfxNormal.xyz[0],
				userGfxNormal.xyz[1], userGfxNormal.xyz[2]);
		break;

	case USER_GFX_VERTEX:
		if (pdbrunOutputVersion < 6)
			count = 0;
		else
			count = sprintf(buf, fmt, userGfxVertex.xyz[0],
				userGfxVertex.xyz[1], userGfxVertex.xyz[2]);
		break;

	case USER_GFX_FONT:
		if (pdbrunOutputVersion < 6)
			count = ::sprintf(buf, fmt, userGfxFont.name,
				userGfxFont.size);
		else
			count = sprintf(buf, fmt, userGfxFont.size,
				userGfxFont.name);
		break;

	case USER_GFX_TEXTPOS:
		if (pdbrunOutputVersion < 6)
			count = 0;
		else
			count = sprintf(buf, fmt, userGfxTextPos.xyz[0],
				userGfxTextPos.xyz[1], userGfxTextPos.xyz[2]);
		break;

	case USER_GFX_LABEL:
		if (pdbrunOutputVersion < 6)
			count = ::sprintf(buf, fmt, userGfxLabel.xyz[0],
				userGfxLabel.xyz[1], userGfxLabel.xyz[2],
				userGfxLabel.text);
		else
			count = sprintf(buf, fmt, userGfxLabel.text);
		break;

	case USER_GFX_MOVE:
		if (pdbrunOutputVersion >= 6)
			count = 0;
		else
			count = sprintf(buf, fmt, userGfxMove.xyz[0],
				userGfxMove.xyz[1], userGfxMove.xyz[2]);
		break;

	case USER_GFX_DRAW:
		if (pdbrunOutputVersion >= 6)
			count = 0;
		else
			count = sprintf(buf, fmt, userGfxDraw.xyz[0],
				userGfxDraw.xyz[1], userGfxDraw.xyz[2]);
		break;

	case USER_GFX_MARKER:
		if (pdbrunOutputVersion >= 6)
			count = 0;
		else
			count = sprintf(buf, fmt, userGfxMarker.xyz[0],
				userGfxMarker.xyz[1], userGfxMarker.xyz[2]);
		break;

	case USER_GFX_POINT:
		if (pdbrunOutputVersion >= 6)
			count = 0;
		else
			count = sprintf(buf, fmt, userGfxPoint.xyz[0],
				userGfxPoint.xyz[1], userGfxPoint.xyz[2]);
		break;

	default:
		count = sprintf(buf, "unknown pdb record #%d", rType);
		break;
	}

	// find last non-blank in buf, and shorten it
	while (count > 1 && isspace(buf[count - 1]))
		count -= 1;
	buf[count] = '\0';
	return buf;
}
