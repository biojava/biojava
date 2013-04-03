// $Id: NHXtags.java,v 1.9 2008/10/05 08:43:22 cmzmasek Exp $
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

package org.forester.io.parsers.nhx;

public final class NHXtags {

    public static final String CUSTOM_DATA_ON_NODE        = "XN=";
    public static final String COLOR                      = "C=";
    public static final String PARENT_BRANCH_WIDTH        = "W=";
    public static final String SUBTREE_NEIGHBORS          = "SNn=";
    public static final String SUPER_ORTHOLOGOUS          = "SOn=";
    public static final String ORTHOLOGOUS                = "On=";
    public static final String TAXONOMY_ID                = "T=";
    public static final String SUPPORT                    = "B=";
    public static final String IS_DUPLICATION             = "D=";
    public static final String ANNOTATION                 = "AN="; //TODO fix on website NHXv2
    public static final String SPECIES_NAME               = "S=";
    public static final String DOMAIN_STRUCTURE           = "DS=";
    public static final String GENE_NAME                  = "GN=";
    public static final String GENE_NAME_SYNONYM          = "G=";
    public static final String SEQUENCE_ACCESSION         = "AC=";
    public static final String NODE_IDENTIFIER            = "ID="; //TODO fix on website NHXv2
    public static final Object BRANCH_WIDTH               = "W=";
    @Deprecated
    public static final String BINARY_DOMAIN_COMBINATIONS = "GDC=";
    @Deprecated
    public static final String DOMAINS_SEPARATOR          = "\\|";
    @Deprecated
    public static final String DOMAINS                    = "GD=";
    @Deprecated
    public static final String EC_NUMBER                  = "E=";
}
