// $Id: PhylogenyDecorator.java,v 1.14 2009/12/09 00:58:22 cmzmasek Exp $
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

package org.forester.tools;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.io.parsers.nhx.NHXFormatException;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Annotation;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.ForesterUtil;

public final class PhylogenyDecorator {

    // From evoruby/lib/evo/apps/tseq_taxonomy_processor.rb:
    final private static String  TP_TAXONOMY_CODE                   = "TAXONOMY_CODE";
    final private static String  TP_TAXONOMY_ID                     = "TAXONOMY_ID";
    final private static String  TP_TAXONOMY_ID_PROVIDER            = "TAXONOMY_ID_PROVIDER";
    final private static String  TP_TAXONOMY_SN                     = "TAXONOMY_SN";
    final private static String  TP_TAXONOMY_CN                     = "TAXONOMY_CN";
    final private static String  TP_TAXONOMY_SYN                    = "TAXONOMY_SYN";
    final private static String  TP_SEQ_SYMBOL                      = "SEQ_SYMBOL";
    final private static String  TP_SEQ_ACCESSION                   = "SEQ_ACCESSION";
    final private static String  TP_SEQ_ACCESSION_SOURCE            = "SEQ_ACCESSION_SOURCE";
    final private static String  TP_SEQ_ANNOTATION_DESC             = "SEQ_ANNOTATION_DESC";
    final private static String  TP_SEQ_ANNOTATION_REF              = "SEQ_ANNOTATION_REF";
    final private static String  TP_SEQ_MOL_SEQ                     = "SEQ_MOL_SEQ";
    final private static String  TP_SEQ_NAME                        = "SEQ_NAME";
    final private static String  TP_NODE_NAME                       = "NODE_NAME";
    final private static Pattern NODENAME_SEQNUMBER_TAXDOMAINNUMBER = Pattern
                                                                            .compile( "^([a-fA-Z0-9]{1,5})_([A-Z0-9]{2,4}[A-Z])(\\d{1,4})$" );
    public final static boolean  SANITIZE                           = false;
    public final static boolean  VERBOSE                            = true;

    private PhylogenyDecorator() {
        // Not needed.
    }

    public static void decorate( final Phylogeny phylogeny,
                                 final Map<String, Map<String, String>> map,
                                 final boolean picky,
                                 final int numbers_of_chars_allowed_to_remove_if_not_found_in_map )
            throws IllegalArgumentException {
        for( final PhylogenyNodeIterator iter = phylogeny.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            final String name = node.getNodeName();
            if ( !ForesterUtil.isEmpty( name ) ) {
                if ( map.containsKey( name ) || ( numbers_of_chars_allowed_to_remove_if_not_found_in_map > 0 ) ) {
                    Map<String, String> new_values = map.get( name );
                    int x = 0;
                    while ( ( new_values == null ) && ( numbers_of_chars_allowed_to_remove_if_not_found_in_map > 0 )
                            && ( x <= numbers_of_chars_allowed_to_remove_if_not_found_in_map ) ) {
                        new_values = map.get( name.substring( 0, name.length() - x ) );
                        ++x;
                    }
                    if ( new_values != null ) {
                        if ( new_values.containsKey( TP_TAXONOMY_CODE ) ) {
                            ForesterUtil.ensurePresenceOfTaxonomy( node );
                            node.getNodeData().getTaxonomy().setTaxonomyCode( new_values.get( TP_TAXONOMY_CODE ) );
                        }
                        if ( new_values.containsKey( TP_TAXONOMY_ID )
                                && new_values.containsKey( TP_TAXONOMY_ID_PROVIDER ) ) {
                            ForesterUtil.ensurePresenceOfTaxonomy( node );
                            node.getNodeData().getTaxonomy().setIdentifier( new Identifier( new_values
                                    .get( TP_TAXONOMY_ID ), new_values.get( TP_TAXONOMY_ID_PROVIDER ) ) );
                        }
                        else if ( new_values.containsKey( TP_TAXONOMY_ID ) ) {
                            ForesterUtil.ensurePresenceOfTaxonomy( node );
                            node.getNodeData().getTaxonomy().setIdentifier( new Identifier( new_values
                                    .get( TP_TAXONOMY_ID ) ) );
                        }
                        if ( new_values.containsKey( TP_TAXONOMY_SN ) ) {
                            ForesterUtil.ensurePresenceOfTaxonomy( node );
                            node.getNodeData().getTaxonomy().setScientificName( new_values.get( TP_TAXONOMY_SN ) );
                        }
                        if ( new_values.containsKey( TP_TAXONOMY_CN ) ) {
                            ForesterUtil.ensurePresenceOfTaxonomy( node );
                            node.getNodeData().getTaxonomy().setCommonName( new_values.get( TP_TAXONOMY_CN ) );
                        }
                        if ( new_values.containsKey( TP_TAXONOMY_SYN ) ) {
                            ForesterUtil.ensurePresenceOfTaxonomy( node );
                            node.getNodeData().getTaxonomy().getSynonyms().add( new_values.get( TP_TAXONOMY_SYN ) );
                        }
                        if ( new_values.containsKey( TP_SEQ_ACCESSION )
                                && new_values.containsKey( TP_SEQ_ACCESSION_SOURCE ) ) {
                            ForesterUtil.ensurePresenceOfSequence( node );
                            node.getNodeData().getSequence().setAccession( new Accession( new_values
                                    .get( TP_SEQ_ACCESSION ), new_values.get( TP_SEQ_ACCESSION_SOURCE ) ) );
                        }
                        if ( new_values.containsKey( TP_SEQ_ANNOTATION_DESC ) ) {
                            ForesterUtil.ensurePresenceOfSequence( node );
                            final Annotation ann = new Annotation();
                            ann.setDesc( new_values.get( TP_SEQ_ANNOTATION_DESC ) );
                            node.getNodeData().getSequence().addAnnotation( ann );
                        }
                        if ( new_values.containsKey( TP_SEQ_ANNOTATION_REF ) ) {
                            ForesterUtil.ensurePresenceOfSequence( node );
                            final Annotation ann = new Annotation();
                            ann.setRef( new_values.get( TP_SEQ_ANNOTATION_REF ) );
                            node.getNodeData().getSequence().addAnnotation( ann );
                        }
                        if ( new_values.containsKey( TP_SEQ_SYMBOL ) ) {
                            ForesterUtil.ensurePresenceOfSequence( node );
                            node.getNodeData().getSequence().setSymbol( new_values.get( TP_SEQ_SYMBOL ) );
                        }
                        if ( new_values.containsKey( TP_SEQ_NAME ) ) {
                            ForesterUtil.ensurePresenceOfSequence( node );
                            node.getNodeData().getSequence().setName( new_values.get( TP_SEQ_NAME ) );
                        }
                        if ( new_values.containsKey( TP_SEQ_MOL_SEQ ) ) {
                            ForesterUtil.ensurePresenceOfSequence( node );
                            node.getNodeData().getSequence().setMolecularSequence( new_values.get( TP_SEQ_MOL_SEQ ) );
                        }
                        if ( new_values.containsKey( TP_NODE_NAME ) ) {
                            node.setName( new_values.get( TP_NODE_NAME ) );
                        }
                    }
                }
                else if ( picky ) {
                    throw new IllegalArgumentException( "\"" + name + "\" not found in name map" );
                }
            }
        }
    }

    /**
     * 
     * 
     * 
     * 
     * 
     * @param phylogeny
     * @param map
     *            maps names (in phylogeny) to new values
     * @param field
     * @param picky
     * @throws IllegalArgumentException
     * @throws NHXFormatException
     */
    public static void decorate( final Phylogeny phylogeny,
                                 final Map<String, String> map,
                                 final FIELD field,
                                 final boolean extract_bracketed_scientific_name,
                                 final boolean picky,
                                 final boolean cut_name_after_space,
                                 final boolean process_name_intelligently,
                                 final boolean process_similar_to,
                                 final int numbers_of_chars_allowed_to_remove_if_not_found_in_map,
                                 final boolean move_domain_numbers_at_end_to_middle ) throws IllegalArgumentException,
            NHXFormatException {
        PhylogenyDecorator.decorate( phylogeny,
                                     map,
                                     field,
                                     extract_bracketed_scientific_name,
                                     picky,
                                     null,
                                     cut_name_after_space,
                                     process_name_intelligently,
                                     process_similar_to,
                                     numbers_of_chars_allowed_to_remove_if_not_found_in_map,
                                     move_domain_numbers_at_end_to_middle );
    }

    /**
     * 
     * 
     * 
     * @param phylogeny
     * @param map
     *            maps names (in phylogeny) to new values if intermediate_map is
     *            null otherwise maps intermediate value to new value
     * @param field
     * @param picky
     * @param intermediate_map
     *            maps name (in phylogeny) to a intermediate value
     * @throws IllegalArgumentException
     */
    public static void decorate( final Phylogeny phylogeny,
                                 final Map<String, String> map,
                                 final FIELD field,
                                 final boolean extract_bracketed_scientific_name,
                                 final boolean picky,
                                 final Map<String, String> intermediate_map,
                                 final boolean cut_name_after_space,
                                 final boolean process_name_intelligently,
                                 final boolean process_similar_to,
                                 final int numbers_of_chars_allowed_to_remove_if_not_found_in_map,
                                 final boolean move_domain_numbers_at_end_to_middle ) throws IllegalArgumentException {
        if ( extract_bracketed_scientific_name && ( field == FIELD.TAXONOMY_SCIENTIFIC_NAME ) ) {
            throw new IllegalArgumentException( "Attempt to extract bracketed scientific name together with data field pointing to scientific name" );
        }
        for( final PhylogenyNodeIterator iter = phylogeny.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            String name = node.getNodeName();
            if ( !ForesterUtil.isEmpty( name ) ) {
                if ( intermediate_map != null ) {
                    name = PhylogenyDecorator.extractIntermediate( intermediate_map, name );
                }
                if ( map.containsKey( name ) || ( numbers_of_chars_allowed_to_remove_if_not_found_in_map > 0 ) ) {
                    String new_value = map.get( name );
                    int x = 0;
                    while ( ( new_value == null ) && ( numbers_of_chars_allowed_to_remove_if_not_found_in_map > 0 )
                            && ( x <= numbers_of_chars_allowed_to_remove_if_not_found_in_map ) ) {
                        new_value = map.get( name.substring( 0, name.length() - x ) );
                        ++x;
                    }
                    if ( new_value != null ) {
                        new_value = new_value.trim();
                        new_value.replaceAll( "/\\s+/", " " );
                        if ( extract_bracketed_scientific_name && new_value.endsWith( "]" ) ) {
                            extractBracketedScientificNames( node, new_value );
                        }
                        switch ( field ) {
                            case SEQUENCE_ANNOTATION_DESC:
                                if ( PhylogenyDecorator.VERBOSE ) {
                                    System.out.println( name + ": " + new_value );
                                }
                                if ( !node.getNodeData().isHasSequence() ) {
                                    node.getNodeData().setSequence( new Sequence() );
                                }
                                final Annotation annotation = new Annotation();
                                annotation.setDesc( new_value );
                                node.getNodeData().getSequence().addAnnotation( annotation );
                                break;
                            case DOMAIN_STRUCTURE:
                                if ( PhylogenyDecorator.VERBOSE ) {
                                    System.out.println( name + ": " + new_value );
                                }
                                if ( !node.getNodeData().isHasSequence() ) {
                                    node.getNodeData().setSequence( new Sequence() );
                                }
                                node.getNodeData().getSequence()
                                        .setDomainArchitecture( new DomainArchitecture( new_value ) );
                                break;
                            case TAXONOMY_CODE:
                                if ( PhylogenyDecorator.VERBOSE ) {
                                    System.out.println( name + ": " + new_value );
                                }
                                ForesterUtil.ensurePresenceOfTaxonomy( node );
                                node.getNodeData().getTaxonomy().setTaxonomyCode( new_value );
                                break;
                            case TAXONOMY_SCIENTIFIC_NAME:
                                if ( PhylogenyDecorator.VERBOSE ) {
                                    System.out.println( name + ": " + new_value );
                                }
                                ForesterUtil.ensurePresenceOfTaxonomy( node );
                                node.getNodeData().getTaxonomy().setScientificName( new_value );
                                break;
                            case SEQUENCE_NAME:
                                if ( PhylogenyDecorator.VERBOSE ) {
                                    System.out.println( name + ": " + new_value );
                                }
                                if ( !node.getNodeData().isHasSequence() ) {
                                    node.getNodeData().setSequence( new Sequence() );
                                }
                                node.getNodeData().getSequence().setName( new_value );
                                break;
                            case NODE_NAME:
                                if ( PhylogenyDecorator.VERBOSE ) {
                                    System.out.print( name + " -> " );
                                }
                                if ( cut_name_after_space ) {
                                    if ( PhylogenyDecorator.VERBOSE ) {
                                        System.out.print( new_value + " -> " );
                                    }
                                    new_value = PhylogenyDecorator.deleteAtFirstSpace( new_value );
                                }
                                else if ( process_name_intelligently ) {
                                    if ( PhylogenyDecorator.VERBOSE ) {
                                        System.out.print( new_value + " -> " );
                                    }
                                    new_value = PhylogenyDecorator.processNameIntelligently( new_value );
                                }
                                else if ( process_similar_to ) {
                                    if ( PhylogenyDecorator.VERBOSE ) {
                                        System.out.print( new_value + " -> " );
                                    }
                                    new_value = PhylogenyDecorator.processSimilarTo( new_value );
                                }
                                if ( PhylogenyDecorator.SANITIZE ) {
                                    new_value = PhylogenyDecorator.sanitize( new_value );
                                }
                                if ( PhylogenyDecorator.VERBOSE ) {
                                    System.out.println( new_value );
                                }
                                node.setName( new_value );
                                break;
                            default:
                                throw new IllegalStateException( "unknown field \"" + field + "\"" );
                        }
                        if ( move_domain_numbers_at_end_to_middle && ( field != FIELD.NODE_NAME ) ) {
                            node.setName( moveDomainNumbersAtEnd( node.getNodeName() ) );
                        }
                    }
                }
                else if ( picky ) {
                    throw new IllegalArgumentException( "\"" + name + "\" not found in name map" );
                }
            }
        }
    }

    public static void decorate( final Phylogeny[] phylogenies,
                                 final Map<String, Map<String, String>> map,
                                 final boolean picky,
                                 final int numbers_of_chars_allowed_to_remove_if_not_found_in_map )
            throws IllegalArgumentException, NHXFormatException {
        for( int i = 0; i < phylogenies.length; ++i ) {
            PhylogenyDecorator.decorate( phylogenies[ i ],
                                         map,
                                         picky,
                                         numbers_of_chars_allowed_to_remove_if_not_found_in_map );
        }
    }

    public static void decorate( final Phylogeny[] phylogenies,
                                 final Map<String, String> map,
                                 final FIELD field,
                                 final boolean extract_bracketed_scientific_name,
                                 final boolean picky,
                                 final boolean cut_name_after_space,
                                 final boolean process_name_intelligently,
                                 final boolean process_similar_to,
                                 final int numbers_of_chars_allowed_to_remove_if_not_found_in_map,
                                 final boolean move_domain_numbers_at_end_to_middle ) throws IllegalArgumentException,
            NHXFormatException {
        for( int i = 0; i < phylogenies.length; ++i ) {
            PhylogenyDecorator.decorate( phylogenies[ i ],
                                         map,
                                         field,
                                         extract_bracketed_scientific_name,
                                         picky,
                                         cut_name_after_space,
                                         process_name_intelligently,
                                         process_similar_to,
                                         numbers_of_chars_allowed_to_remove_if_not_found_in_map,
                                         move_domain_numbers_at_end_to_middle );
        }
    }

    public static void decorate( final Phylogeny[] phylogenies,
                                 final Map<String, String> map,
                                 final FIELD field,
                                 final boolean extract_bracketed_scientific_name,
                                 final boolean picky,
                                 final Map<String, String> intermediate_map,
                                 final boolean cut_name_after_space,
                                 final boolean process_name_intelligently,
                                 final boolean process_similar_to,
                                 final int numbers_of_chars_allowed_to_remove_if_not_found_in_map,
                                 final boolean move_domain_numbers_at_end_to_middle ) throws IllegalArgumentException,
            NHXFormatException {
        for( int i = 0; i < phylogenies.length; ++i ) {
            PhylogenyDecorator.decorate( phylogenies[ i ],
                                         map,
                                         field,
                                         extract_bracketed_scientific_name,
                                         picky,
                                         intermediate_map,
                                         cut_name_after_space,
                                         process_name_intelligently,
                                         process_similar_to,
                                         numbers_of_chars_allowed_to_remove_if_not_found_in_map,
                                         move_domain_numbers_at_end_to_middle );
        }
    }

    private static String deleteAtFirstSpace( final String name ) {
        final int first_space = name.indexOf( " " );
        if ( first_space > 1 ) {
            return name.substring( 0, first_space ).trim();
        }
        return name;
    }

    private static void extractBracketedScientificNames( final PhylogenyNode node, final String new_value ) {
        final int i = new_value.lastIndexOf( "[" );
        final String scientific_name = new_value.substring( i + 1, new_value.length() - 1 );
        ForesterUtil.ensurePresenceOfTaxonomy( node );
        node.getNodeData().getTaxonomy().setScientificName( scientific_name );
    }

    private static String extractIntermediate( final Map<String, String> intermediate_map, final String name ) {
        String new_name = null;
        if ( PhylogenyDecorator.VERBOSE ) {
            System.out.print( name + " => " );
        }
        if ( intermediate_map.containsKey( name ) ) {
            new_name = intermediate_map.get( name );
            if ( ForesterUtil.isEmpty( new_name ) ) {
                throw new IllegalArgumentException( "\"" + name + "\" maps to null or empty string in secondary map" );
            }
        }
        else {
            throw new IllegalArgumentException( "\"" + name + "\" not found in name secondary map" );
        }
        if ( PhylogenyDecorator.VERBOSE ) {
            System.out.println( new_name + "  " );
        }
        return new_name;
    }

    private static String moveDomainNumbersAtEnd( final String node_name ) {
        final Matcher m = NODENAME_SEQNUMBER_TAXDOMAINNUMBER.matcher( node_name );
        if ( m.matches() ) {
            final String seq_number = m.group( 1 );
            final String tax = m.group( 2 );
            final String domain_number = m.group( 3 );
            return seq_number + "_[" + domain_number + "]_" + tax;
        }
        else {
            return node_name;
        }
    }

    public static Map<String, Map<String, String>> parseMappingTable( final File mapping_table_file )
            throws IOException {
        final Map<String, Map<String, String>> map = new HashMap<String, Map<String, String>>();
        BasicTable<String> mapping_table = null;
        mapping_table = BasicTableParser.parse( mapping_table_file, "\t", false );
        for( int row = 0; row < mapping_table.getNumberOfRows(); ++row ) {
            final Map<String, String> row_map = new HashMap<String, String>();
            String name = null;
            for( int col = 0; col < mapping_table.getNumberOfColumns(); ++col ) {
                final String table_cell = mapping_table.getValue( col, row );
                if ( col == 0 ) {
                    name = table_cell;
                }
                else if ( table_cell != null ) {
                    final String key = table_cell.substring( 0, table_cell.indexOf( ':' ) );
                    final String val = table_cell.substring( table_cell.indexOf( ':' ) + 1, table_cell.length() );
                    row_map.put( key, val );
                }
            }
            map.put( name, row_map );
        }
        return map;
    }

    private static String processNameIntelligently( final String name ) {
        final String[] s = name.split( " " );
        if ( s.length < 2 ) {
            return name;
        }
        else if ( ( s[ 0 ].indexOf( "_" ) > 0 ) && ( s[ 0 ].indexOf( "|" ) > 0 ) ) {
            return s[ 0 ];
        }
        else if ( ( s[ 1 ].indexOf( "_" ) > 0 ) && ( s[ 1 ].indexOf( "|" ) > 0 ) ) {
            return s[ 1 ];
        }
        else if ( ( s[ 0 ].indexOf( "_" ) > 0 ) && ( s[ 0 ].indexOf( "." ) > 0 ) ) {
            return s[ 0 ];
        }
        else if ( ( s[ 1 ].indexOf( "_" ) > 0 ) && ( s[ 1 ].indexOf( "." ) > 0 ) ) {
            return s[ 1 ];
        }
        else if ( s[ 0 ].indexOf( "_" ) > 0 ) {
            return s[ 0 ];
        }
        else if ( s[ 1 ].indexOf( "_" ) > 0 ) {
            return s[ 1 ];
        }
        else {
            return s[ 0 ];
        }
    }

    private static String processSimilarTo( final String name ) {
        final int i = name.toLowerCase().indexOf( "similar to" );
        String similar_to = "";
        if ( i >= 0 ) {
            similar_to = " similarity=" + name.substring( i + 10 ).trim();
        }
        final String pi = processNameIntelligently( name );
        return pi + similar_to;
    }

    private static String sanitize( String s ) {
        s = s.replace( ' ', '_' );
        s = s.replace( '(', '{' );
        s = s.replace( ')', '}' );
        s = s.replace( '[', '{' );
        s = s.replace( ']', '}' );
        s = s.replace( ',', '_' );
        return s;
    }

    public static enum FIELD {
        NODE_NAME, SEQUENCE_ANNOTATION_DESC, DOMAIN_STRUCTURE, TAXONOMY_CODE, TAXONOMY_SCIENTIFIC_NAME, SEQUENCE_NAME;
    }
}
