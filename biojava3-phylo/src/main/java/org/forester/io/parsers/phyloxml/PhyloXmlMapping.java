// $Id: PhyloXmlMapping.java,v 1.9 2009/11/20 22:22:10 cmzmasek Exp $
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

package org.forester.io.parsers.phyloxml;

/*
 * @author Christian Zmasek TODO To change the template for this generated type
 * comment go to Window - Preferences - Java - Code Style - Code Templates
 */
public final class PhyloXmlMapping {

    public static final String PHYLOGENY                                           = "phylogeny";
    public static final String PHYLOGENY_NAME                                      = "name";
    public static final String PHYLOGENY_DESCRIPTION                               = "description";
    public static final String PHYLOGENY_IS_REROOTABLE_ATTR                        = "rerootable";
    public static final String PHYLOGENY_BRANCHLENGTH_UNIT_ATTR                    = "branch_length_unit";
    public static final String PHYLOGENY_IS_ROOTED_ATTR                            = "rooted";
    public static final String PHYLOGENY_TYPE_ATTR                                 = "type";
    public static final String CLADE                                               = "clade";
    public static final String NODE_NAME                                           = "name";
    public static final String SEQUENCE                                            = "sequence";
    public static final String SEQUENCE_NAME                                       = "name";
    public static final String SEQUENCE_SYMBOL                                     = "symbol";
    public static final String ACCESSION                                           = "accession";
    public static final String ACCESSION_SOURCE_ATTR                               = "source";
    public static final String SEQUENCE_LOCATION                                   = "location";
    public static final String SEQUENCE_MOL_SEQ                                    = "mol_seq";
    public static final String ANNOTATION                                          = "annotation";
    public static final String ANNOTATION_DESC                                     = "desc";
    public static final String ANNOTATION_REF_ATTR                                 = "ref";
    public static final String ANNOTATION_EVIDENCE_ATTR                            = "evidence";
    public static final String ANNOTATION_TYPE_ATTR                                = "type";
    public static final String TAXONOMY                                            = "taxonomy";
    public static final String TAXONOMY_SCIENTIFIC_NAME                            = "scientific_name";
    public static final String TAXONOMY_COMMON_NAME                                = "common_name";
    public static final String TAXONOMY_CODE                                       = "code";
    public static final String TAXONOMY_RANK                                       = "rank";
    public static final String TAXONOMY_SYNONYM                                    = "synonym";
    public static final String TAXONOMY_AUTHORITY                                  = "authority";
    public static final String DISTRIBUTION                                        = "distribution";
    public static final String BINARY_CHARACTERS                                   = "binary_characters";
    public static final String BINARY_CHARACTERS_PRESENT                           = "present";
    public static final String BINARY_CHARACTERS_GAINED                            = "gained";
    public static final String BINARY_CHARACTERS_LOST                              = "lost";
    public static final String BINARY_CHARACTERS_TYPE_ATTR                         = "type";
    public static final String BINARY_CHARACTERS_PRESENT_COUNT_ATTR                = "present_count";
    public static final String BINARY_CHARACTERS_GAINED_COUNT_ATTR                 = "gained_count";
    public static final String BINARY_CHARACTERS_LOST_COUNT_ATTR                   = "lost_count";
    public static final String BRANCH_LENGTH                                       = "branch_length";
    public static final String CONFIDENCE                                          = "confidence";
    public static final String CONFIDENCE_TYPE_ATTR                                = "type";
    public static final String COLOR                                               = "color";
    public static final String COLOR_RED                                           = "red";
    public static final String COLOR_GREEN                                         = "green";
    public static final String COLOR_BLUE                                          = "blue";
    public final static String SEQUENCE_DOMAIN_ARCHITECTURE_DOMAIN                 = "domain";
    public final static String SEQUENCE_DOMAIN_ARCHITECTURE_PROT_DOMAIN_FROM       = "from";
    public final static String SEQUENCE_DOMAIN_ARCHITECTURE_PROT_DOMAIN_TO         = "to";
    public final static String SEQUENCE_DOMAIN_ARCHITECTURE_PROT_DOMAIN_CONFIDENCE = "confidence";
    // public final static String NODE_IDENTIFIER                                     = "node_id";
    public final static String IDENTIFIER                                          = "id";
    public final static String IDENTIFIER_PROVIDER_ATTR                            = "provider";
    public static final String URI                                                 = "uri";
    public static final String WIDTH                                               = "width";
    public final static String EVENTS                                              = "events";
    public final static String EVENT_TYPE                                          = "type";
    public final static String EVENT_DUPLICATIONS                                  = "duplications";
    public final static String EVENT_SPECIATIONS                                   = "speciations";
    public final static String EVENT_LOSSES                                        = "losses";
    public final static String SEQUENCE_DOMAIN_ARCHITECURE                         = "domain_architecture";
    public final static String SEQUENCE_DOMAIN_ARCHITECTURE_LENGTH                 = "length";
    public final static String SEQUENCE_TYPE                                       = "type";
    public static final String BINARY_CHARACTER                                    = "bc";
    public static final String URI_DESC_ATTR                                       = "desc";
    public static final String TYPE_ATTR                                           = "type";
    public static final String REFERENCE                                           = "reference";
    public static final String REFERENCE_DOI_ATTR                                  = "doi";
    public static final String REFERENCE_DESC                                      = "desc";
    public static final String PROPERTY                                            = "property";
    public static final String PROPERTY_REF                                        = "ref";
    public static final String PROPERTY_UNIT                                       = "unit";
    public static final String PROPERTY_DATATYPE                                   = "datatype";
    public static final String PROPERTY_APPLIES_TO                                 = "applies_to";
    public static final String ID_REF                                              = "id_ref";
    public static final String ANNOTATION_SOURCE_ATTR                              = "source";
    public static final String DISTRIBUTION_DESC                                   = "desc";
    public static final String POINT                                               = "point";
    public static final String POINT_LONGITUDE                                     = "longitude";
    public static final String POINT_LATITUDE                                      = "latitude";
    public static final String POINT_ALTITUDE                                      = "alt";
    public static final String POINT_ALTITUDE_UNIT_ATTR                            = "alt_unit";
    public static final String POINT_GEODETIC_DATUM                                = "geodetic_datum";
    public static final String CLADE_DATE                                          = "date";
    public static final String CLADE_DATE_UNIT                                     = "unit";
    public static final String CLADE_DATE_DESC                                     = "desc";
    public static final String CLADE_DATE_MIN                                      = "minimum";
    public static final String CLADE_DATE_MAX                                      = "maximum";
    public static final String CLADE_DATE_VALUE                                    = "value";

    /**
     * 
     */
    private PhyloXmlMapping() {
    }
}
