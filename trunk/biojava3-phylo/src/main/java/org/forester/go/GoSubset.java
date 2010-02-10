// $Id: GoSubset.java,v 1.3 2009/11/10 19:57:09 cmzmasek Exp $
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

public interface GoSubset extends Comparable<GoSubset> {

    public static final String GOSLIM_GENERIC_STR = "goslim_generic";
    public static final String GOSLIM_GOA_STR     = "goslim_goa";
    public static final String GOSLIM_PIR_STR     = "goslim_pir";
    public static final String GOSUBSET_PROK_STR  = "gosubset_prok";
    public static final String GOSLIM_CANDIDA_STR = "goslim_candida";
    public static final String GOSLIM_PLANT_STR   = "goslim_plant";
    public static final String GOSLIM_YEAST_STR   = "goslim_yeast";
    public static final String GOSLIM_POMBE_STR   = "goslim_pombe";

    public Type getType();

    public static enum Type {
        GOSLIM_GENERIC, GOSLIM_GOA, GOSLIM_PIR, GOSUBSET_PROK, GOSLIM_CANDIDA, GOSLIM_PLANT, GOSLIM_YEAST, GOSLIM_POMBE;
    }
}
