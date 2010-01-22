// $Id: PhyloXmlUtil.java,v 1.2 2009/12/17 02:28:13 cmzmasek Exp $
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

package org.forester.io.parsers.phyloxml;

import java.util.HashSet;
import java.util.Set;
import java.util.regex.Pattern;

public final class PhyloXmlUtil {

    public final static Pattern     SEQUENCE_SYMBOL_PATTERN                    = Pattern.compile( "\\S{1,20}" );
    public final static Pattern     TAXOMONY_CODE_PATTERN                      = Pattern.compile( "[a-zA-Z0-9_]{1,10}" );
    public final static Pattern     LIT_REF_DOI_PATTERN                        = Pattern
                                                                                       .compile( "[a-zA-Z0-9_\\.]+\\S+" );
    public final static Set<String> SEQUENCE_TYPES                             = new HashSet<String>();
    public final static Set<String> TAXONOMY_RANKS                             = new HashSet<String>();
    public static final int         ROUNDING_DIGITS_FOR_PHYLOXML_DOUBLE_OUTPUT = 9;
    static {
        SEQUENCE_TYPES.add( "rna" );
        SEQUENCE_TYPES.add( "protein" );
        SEQUENCE_TYPES.add( "dna" );
        TAXONOMY_RANKS.add( "domain" );
        TAXONOMY_RANKS.add( "kingdom" );
        TAXONOMY_RANKS.add( "subkingdom" );
        TAXONOMY_RANKS.add( "branch" );
        TAXONOMY_RANKS.add( "infrakingdom" );
        TAXONOMY_RANKS.add( "superphylum" );
        TAXONOMY_RANKS.add( "phylum" );
        TAXONOMY_RANKS.add( "subphylum" );
        TAXONOMY_RANKS.add( "infraphylum" );
        TAXONOMY_RANKS.add( "microphylum" );
        TAXONOMY_RANKS.add( "superdivision" );
        TAXONOMY_RANKS.add( "division" );
        TAXONOMY_RANKS.add( "subdivision" );
        TAXONOMY_RANKS.add( "infradivision" );
        TAXONOMY_RANKS.add( "superclass" );
        TAXONOMY_RANKS.add( "class" );
        TAXONOMY_RANKS.add( "subclass" );
        TAXONOMY_RANKS.add( "infraclass" );
        TAXONOMY_RANKS.add( "superlegion" );
        TAXONOMY_RANKS.add( "legion" );
        TAXONOMY_RANKS.add( "sublegion" );
        TAXONOMY_RANKS.add( "infralegion" );
        TAXONOMY_RANKS.add( "supercohort" );
        TAXONOMY_RANKS.add( "cohort" );
        TAXONOMY_RANKS.add( "subcohort" );
        TAXONOMY_RANKS.add( "infracohort" );
        TAXONOMY_RANKS.add( "superorder" );
        TAXONOMY_RANKS.add( "order" );
        TAXONOMY_RANKS.add( "suborder" );
        TAXONOMY_RANKS.add( "superfamily" );
        TAXONOMY_RANKS.add( "family" );
        TAXONOMY_RANKS.add( "subfamily" );
        TAXONOMY_RANKS.add( "supertribe" );
        TAXONOMY_RANKS.add( "tribe" );
        TAXONOMY_RANKS.add( "subtribe" );
        TAXONOMY_RANKS.add( "infratribe" );
        TAXONOMY_RANKS.add( "genus" );
        TAXONOMY_RANKS.add( "subgenus" );
        TAXONOMY_RANKS.add( "superspecies" );
        TAXONOMY_RANKS.add( "species" );
        TAXONOMY_RANKS.add( "subspecies" );
        TAXONOMY_RANKS.add( "variety" );
        TAXONOMY_RANKS.add( "subvariety" );
        TAXONOMY_RANKS.add( "form" );
        TAXONOMY_RANKS.add( "subform" );
        TAXONOMY_RANKS.add( "cultivar" );
        TAXONOMY_RANKS.add( "unknown" );
        TAXONOMY_RANKS.add( "other" );
    };
}
