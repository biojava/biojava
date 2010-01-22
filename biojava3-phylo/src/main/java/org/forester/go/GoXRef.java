// $Id: GoXRef.java,v 1.14 2009/11/10 19:57:09 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// All rights reserved
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: cmzmasek@yahoo.com
// WWW: www.phylosoft.org/forester

package org.forester.go;

public interface GoXRef extends Comparable<GoXRef> {

    public static final String EC_STR                = "EC";
    public static final String META_CYC_STR          = "MetaCyc";
    public static final String REACTOME_STR          = "Reactome";
    public static final String RESID_STR             = "RESID";
    public static final String UM_BBD_ENZYME_ID_STR  = "UM-BBD_enzymeID";
    public static final String UM_BBD_PATHWAY_ID_STR = "UM-BBD_pathwayID";
    public static final String UM_BBD_REACTIONID_STR = "UM-BBD_reactionID";
    public static final String TC_STR                = "TC";
    public static final String ARACYC_STR            = "AraCyc";
    public static final String XX_STR                = "XX";
    public static final String PMID_STR              = "PMID";
    public static final String IMG_STR               = "IMG";
    public static final String GOC_STR               = "GOC";
    public static final String WIKIPEDIA_STR         = "Wikipedia";
    public static final String KEGG_STR              = "KEGG";

    public Type getType();

    public String getXRef();

    public static enum Type {
        EC,
        META_CYC,
        REACTOME,
        RESID,
        UM_BBD_ENZYME_ID,
        UM_BBD_PATHWAY_ID,
        UM_BBD_REACTIONID,
        TC,
        ARACYC,
        XX,
        PMID,
        IMG,
        GOC,
        WIKIPEDIA,
        KEGG;
    }
}
