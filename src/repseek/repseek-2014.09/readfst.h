/*
	Copyright (C) 2006 G achaz, F boyer, E rocha, A viari and E coissac.

	This program is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public License
	as published by the Free Software Foundation; either version 2.1
	of the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 
 	for more information, please contact guillaume achaz <achaz@abi.snv.jussieu.fr>
*/
/**
 * @file   KMRK_readfst.h
 * @author Eric Coissac <coissac@inrialpes.fr>
 * @date   Tue Feb 24 22:11:07 2004
 * 
 * @brief  Header file for fasta file reading
 * 
 * 
 */


#ifndef KMRK_readfst_h
#define KMRK_readfst_h

char *file2mem(char *file, int32_t *size);  

char *readFastaSeq(char *file, int32_t *size, int32_t *nX);

char *fuse2files_mem(char *file1, char *file2, int32_t *size1, int32_t *size2);

#endif /* KMRK_readfst_h */
