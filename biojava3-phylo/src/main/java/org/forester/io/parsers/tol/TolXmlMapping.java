// $Id: TolXmlMapping.java,v 1.7 2009/11/17 19:53:23 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// Copyright (C) 2000-2001 Washington University School of Medicine
// and Howard Hughes Medical Institute
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

package org.forester.io.parsers.tol;

public final class TolXmlMapping {

    public static final String PHYLOGENY               = "TREE";
    public static final String CLADE                   = "NODE";
    public static final String AUTHDATE                = "AUTHDATE";
    public static final String AUTHORITY               = "AUTHORITY";
    public static final String TAXONOMY_NAME           = "NAME";
    public static final String OTHERNAMES              = "OTHERNAMES";
    public static final String OTHERNAME               = "OTHERNAME";
    public static final String OTHERNAME_NAME          = "NAME";
    public static final String NODE_ID_ATTR            = "ID";
    public static final String NODE_ITALICIZENAME_ATTR = "ITALICIZENAME";
    public static final String TOL_TAXONOMY_ID_TYPE    = "tol";

    private TolXmlMapping() {
        // Hidden.
    }
}