// $Id: decorator.java,v 1.32 2009/11/20 22:22:09 cmzmasek Exp $
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

package org.forester.application;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.tools.PhylogenyDecorator;
import org.forester.tools.PhylogenyDecorator.FIELD;
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public final class decorator {

    private static final String  SEQUENCE_NAME_FIELD                    = "s";
    private static final String  TAXONOMY_CODE_FIELD                    = "c";
    private static final String  TAXONOMY_SCIENTIFIC_NAME_FIELD         = "sn";
    private static final String  DS_FILED                               = "d";
    private static final String  SEQUENCE_ANNOTATION_DESC               = "a";
    private static final String  NODE_NAME_FIELD                        = "n";
    final static private String  PICKY_OPTION                           = "p";
    final static private String  FIELD_OPTION                           = "f";
    final static private String  MOVE_DOMAIN_NUMBER_OPTION              = "mdn";       // Hidden expert option.
    final static private String  TREE_NAME_OPTION                       = "pn";
    final static private String  TREE_ID_OPTION                         = "pi";
    final static private String  TREE_DESC_OPTION                       = "pd";
    final static private String  EXTRACT_BRACKETED_SCIENTIC_NAME_OPTION = "sn";
    final static private String  PROCESS_NAME_INTELLIGENTLY_OPTION      = "x";
    final static private String  PROCESS_SIMILAR_TO_OPTION              = "xs";
    final static private String  CUT_NAME_AFTER_FIRST_SPACE_OPTION      = "c";
    final static private String  ALLOW_REMOVAL_OF_CHARS_OPTION          = "r";
    final static private String  ADVANCED_TABLE_OPTION                  = "table";
    final static private String  KEY_COLUMN                             = "k";
    final static private String  VALUE_COLUMN                           = "v";
    final static private String  MAPPING_FILE_SEPARATOR_OPTION          = "s";
    final static private String  MAPPING_FILE_SEPARATOR_DEFAULT         = ":";
    final static private boolean USE_FIRST_SEPARATOR_ONLY               = true;
    final static private String  PRG_NAME                               = "decorator";
    final static private String  PRG_VERSION                            = "1.10";
    final static private String  PRG_DATE                               = "2009.10.08";

    private static void argumentsError() {
        System.out.println();
        System.out.println( decorator.PRG_NAME + " -" + ADVANCED_TABLE_OPTION + " | -f=<c> <phylogenies infile> "
                + "<mapping table file> <phylogenies outfile>" );
        System.out.println();
        System.out.println( "options:" );
        System.out.println();
        System.out.println( " -" + ADVANCED_TABLE_OPTION + " : table instead of one to one map (-f=<c>)" );
        System.out.println( " -r=<n> : allow to remove up to n characters from the end of the names" );
        System.out.println( "          in phylogenies infile if not found (in map) otherwise" );
        System.out.println( " -p     : for picky, fails if node name not found in mapping table, default is off" );
        System.out.println( " -" + TREE_NAME_OPTION + "=<s>: name for the phylogeny" );
        System.out.println( " -" + TREE_ID_OPTION + "=<s>: identifier for the phylogeny (in the form provider:value)" );
        System.out.println( " -" + TREE_DESC_OPTION + "=<s>: description for phylogenies" );
        System.out.println();
        System.out.println();
        System.out.println( "advanced options, only available if -" + ADVANCED_TABLE_OPTION + " is not used:" );
        System.out.println();
        System.out.println( " -f=<c> : field to be replaced: " + NODE_NAME_FIELD + " : node name" );
        System.out.println( "                                " + SEQUENCE_ANNOTATION_DESC
                + " : sequence annotation description" );
        System.out.println( "                                " + DS_FILED + " : domain structure" );
        System.out.println( "                                " + TAXONOMY_CODE_FIELD + " : taxonomy code" );
        System.out.println( "                                " + TAXONOMY_SCIENTIFIC_NAME_FIELD
                + ": taxonomy scientific name" );
        System.out.println( "                                " + SEQUENCE_NAME_FIELD + " : sequence name" );
        System.out.println( " -k=<n> : key column in mapping table (0 based)," );
        System.out.println( "          names of the node to be decorated - default is 0" );
        System.out.println( " -v=<n> : value column in mapping table (0 based)," );
        System.out.println( "          data which with to decorate - default is 1" );
        System.out.println( " -" + EXTRACT_BRACKETED_SCIENTIC_NAME_OPTION
                + "    : to extract bracketed scientific names" );
        System.out.println( " -s=<c> : column separator in mapping file, default is \""
                + decorator.MAPPING_FILE_SEPARATOR_DEFAULT + "\"" );
        System.out.println( " -x     : process name \"intelligently\" (only for -f=n)" );
        System.out.println( " -" + decorator.PROCESS_SIMILAR_TO_OPTION
                + "    : process name \"intelligently\" and process information after \"similar to\" (only for -f=n)" );
        System.out.println( " -c     : cut name after first space (only for -f=n)" );
        System.out.println();
        System.exit( -1 );
    }

    public static void main( final String args[] ) {
        ForesterUtil.printProgramInformation( decorator.PRG_NAME, decorator.PRG_VERSION, decorator.PRG_DATE );
        if ( ( args.length < 4 ) || ( args.length > 12 ) ) {
            decorator.argumentsError();
        }
        CommandLineArguments cla = null;
        try {
            cla = new CommandLineArguments( args );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        if ( ( cla.getNumberOfNames() < 3 ) || ( cla.getNumberOfNames() > 4 ) ) {
            decorator.argumentsError();
        }
        final File phylogenies_infile = cla.getFile( 0 );
        final File mapping_infile = cla.getFile( 1 );
        final File phylogenies_outfile = cla.getFile( 2 );
        if ( phylogenies_outfile.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + phylogenies_outfile + "] already exists" );
        }
        final List<String> allowed_options = new ArrayList<String>();
        allowed_options.add( decorator.ADVANCED_TABLE_OPTION );
        allowed_options.add( decorator.PICKY_OPTION );
        allowed_options.add( decorator.FIELD_OPTION );
        allowed_options.add( decorator.PROCESS_NAME_INTELLIGENTLY_OPTION );
        allowed_options.add( decorator.PROCESS_SIMILAR_TO_OPTION );
        allowed_options.add( decorator.CUT_NAME_AFTER_FIRST_SPACE_OPTION );
        allowed_options.add( decorator.ALLOW_REMOVAL_OF_CHARS_OPTION );
        allowed_options.add( decorator.KEY_COLUMN );
        allowed_options.add( decorator.VALUE_COLUMN );
        allowed_options.add( decorator.MAPPING_FILE_SEPARATOR_OPTION );
        allowed_options.add( decorator.EXTRACT_BRACKETED_SCIENTIC_NAME_OPTION );
        allowed_options.add( decorator.TREE_NAME_OPTION );
        allowed_options.add( decorator.TREE_ID_OPTION );
        allowed_options.add( decorator.TREE_DESC_OPTION );
        allowed_options.add( decorator.MOVE_DOMAIN_NUMBER_OPTION );
        final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
        if ( dissallowed_options.length() > 0 ) {
            ForesterUtil.fatalError( decorator.PRG_NAME, "unknown option(s): " + dissallowed_options );
        }
        final boolean advanced_table = cla.isOptionSet( decorator.ADVANCED_TABLE_OPTION );
        if ( !advanced_table ) {
            final List<String> mandatory_options = new ArrayList<String>();
            mandatory_options.add( decorator.FIELD_OPTION );
            final String missing_options = cla.validateMandatoryOptionsAsString( mandatory_options );
            if ( missing_options.length() > 0 ) {
                ForesterUtil.fatalError( decorator.PRG_NAME, "missing option(s): " + missing_options );
            }
        }
        final boolean picky = cla.isOptionSet( decorator.PICKY_OPTION );
        String separator = decorator.MAPPING_FILE_SEPARATOR_DEFAULT;
        if ( cla.isOptionSet( decorator.MAPPING_FILE_SEPARATOR_OPTION ) ) {
            if ( advanced_table ) {
                argumentsError();
            }
            separator = cla.getOptionValue( decorator.MAPPING_FILE_SEPARATOR_OPTION );
        }
        int key_column = 0;
        int value_column = 1;
        String field_str = "";
        FIELD field = FIELD.NODE_NAME;
        int numbers_of_chars_allowed_to_remove_if_not_found_in_map = -1;
        boolean cut_name_after_space = false;
        boolean process_name_intelligently = false;
        boolean process_similar_to = false;
        boolean extract_bracketed_scientific_name = false;
        boolean move_domain_numbers_at_end_to_middle = false;
        String tree_name = "";
        String tree_id = "";
        String tree_desc = "";
        try {
            if ( cla.isOptionSet( decorator.TREE_NAME_OPTION ) ) {
                tree_name = cla.getOptionValueAsCleanString( decorator.TREE_NAME_OPTION );
            }
            if ( cla.isOptionSet( decorator.TREE_ID_OPTION ) ) {
                tree_id = cla.getOptionValueAsCleanString( decorator.TREE_ID_OPTION );
            }
            if ( cla.isOptionSet( decorator.TREE_DESC_OPTION ) ) {
                tree_desc = cla.getOptionValueAsCleanString( decorator.TREE_DESC_OPTION );
            }
            if ( cla.isOptionSet( decorator.EXTRACT_BRACKETED_SCIENTIC_NAME_OPTION ) ) {
                if ( advanced_table ) {
                    argumentsError();
                }
                extract_bracketed_scientific_name = true;
            }
            if ( cla.isOptionSet( decorator.KEY_COLUMN ) ) {
                if ( advanced_table ) {
                    argumentsError();
                }
                key_column = cla.getOptionValueAsInt( decorator.KEY_COLUMN );
            }
            if ( cla.isOptionSet( decorator.VALUE_COLUMN ) ) {
                if ( advanced_table ) {
                    argumentsError();
                }
                value_column = cla.getOptionValueAsInt( decorator.VALUE_COLUMN );
            }
            if ( cla.isOptionSet( decorator.CUT_NAME_AFTER_FIRST_SPACE_OPTION ) ) {
                if ( advanced_table ) {
                    argumentsError();
                }
                cut_name_after_space = true;
            }
            if ( cla.isOptionSet( decorator.PROCESS_NAME_INTELLIGENTLY_OPTION ) ) {
                if ( advanced_table ) {
                    argumentsError();
                }
                process_name_intelligently = true;
            }
            if ( cla.isOptionSet( decorator.PROCESS_SIMILAR_TO_OPTION ) ) {
                if ( advanced_table ) {
                    argumentsError();
                }
                process_similar_to = true;
            }
            if ( cla.isOptionSet( decorator.ALLOW_REMOVAL_OF_CHARS_OPTION ) ) {
                numbers_of_chars_allowed_to_remove_if_not_found_in_map = cla
                        .getOptionValueAsInt( decorator.ALLOW_REMOVAL_OF_CHARS_OPTION );
            }
            if ( cla.isOptionSet( decorator.MOVE_DOMAIN_NUMBER_OPTION ) ) {
                move_domain_numbers_at_end_to_middle = true;
            }
            if ( cla.isOptionSet( decorator.FIELD_OPTION ) ) {
                field_str = cla.getOptionValue( decorator.FIELD_OPTION );
                if ( field_str.equals( NODE_NAME_FIELD ) ) {
                    field = FIELD.NODE_NAME;
                }
                else if ( field_str.equals( SEQUENCE_ANNOTATION_DESC ) ) {
                    field = FIELD.SEQUENCE_ANNOTATION_DESC;
                }
                else if ( field_str.equals( DS_FILED ) ) {
                    field = FIELD.DOMAIN_STRUCTURE;
                    extract_bracketed_scientific_name = false;
                }
                else if ( field_str.equals( TAXONOMY_CODE_FIELD ) ) {
                    field = FIELD.TAXONOMY_CODE;
                }
                else if ( field_str.equals( SEQUENCE_NAME_FIELD ) ) {
                    field = FIELD.SEQUENCE_NAME;
                }
                else if ( field_str.equals( TAXONOMY_SCIENTIFIC_NAME_FIELD ) ) {
                    field = FIELD.TAXONOMY_SCIENTIFIC_NAME;
                    extract_bracketed_scientific_name = false;
                }
                else {
                    ForesterUtil.fatalError( decorator.PRG_NAME, "unknown value for \"" + decorator.FIELD_OPTION
                            + "\" option: \"" + field_str + "\"" );
                }
            }
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( decorator.PRG_NAME, "error in command line: " + e.getMessage() );
        }
        if ( ( field != FIELD.NODE_NAME ) && ( cut_name_after_space || process_name_intelligently ) ) {
            ForesterUtil.fatalError( decorator.PRG_NAME, "attempt to use -x or -c option without -f=n" );
        }
        if ( ( field != FIELD.NODE_NAME ) && process_similar_to ) {
            ForesterUtil.fatalError( decorator.PRG_NAME, "attempt to use -" + decorator.PROCESS_SIMILAR_TO_OPTION
                    + " option without -f=n" );
        }
        if ( cut_name_after_space && process_name_intelligently ) {
            ForesterUtil.fatalError( decorator.PRG_NAME, "attempt to use -x and -c option together" );
        }
        if ( process_similar_to && process_name_intelligently ) {
            ForesterUtil.fatalError( decorator.PRG_NAME, "attempt to use -" + decorator.PROCESS_SIMILAR_TO_OPTION
                    + " and -x option together" );
        }
        if ( process_similar_to && cut_name_after_space ) {
            ForesterUtil.fatalError( decorator.PRG_NAME, "attempt to use -" + decorator.PROCESS_SIMILAR_TO_OPTION
                    + " and -c option together" );
        }
        Phylogeny[] phylogenies = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ForesterUtil.createParserDependingOnFileType( phylogenies_infile, true );
            phylogenies = factory.create( phylogenies_infile, pp );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( decorator.PRG_NAME, "failed to read phylgenies from [" + phylogenies_infile
                    + "] [" + e.getMessage() + "]" );
        }
        Map<String, String> map = null;
        if ( !advanced_table ) {
            BasicTable<String> mapping_table = null;
            try {
                mapping_table = BasicTableParser.parse( mapping_infile, separator, decorator.USE_FIRST_SEPARATOR_ONLY );
            }
            catch ( final Exception e ) {
                ForesterUtil.fatalError( decorator.PRG_NAME, "failed to read [" + mapping_infile + "] ["
                        + e.getMessage() + "]" );
            }
            if ( ( key_column < 0 ) || ( key_column >= mapping_table.getNumberOfColumns() ) ) {
                ForesterUtil.fatalError( decorator.PRG_NAME, "illegal value for key column" );
            }
            if ( ( value_column < 0 ) || ( value_column >= mapping_table.getNumberOfColumns() ) ) {
                ForesterUtil.fatalError( decorator.PRG_NAME, "illegal value for value column" );
            }
            map = mapping_table.getColumnsAsMap( key_column, value_column );
        }
        if ( !ForesterUtil.isEmpty( tree_name ) || !ForesterUtil.isEmpty( tree_id )
                || !ForesterUtil.isEmpty( tree_desc ) ) {
            if ( ( phylogenies.length > 1 )
                    && ( !ForesterUtil.isEmpty( tree_name ) || !ForesterUtil.isEmpty( tree_id ) ) ) {
                ForesterUtil.fatalError( decorator.PRG_NAME,
                                         "attempt to set same name or id on more than one phylogeny" );
            }
            if ( !ForesterUtil.isEmpty( tree_name ) ) {
                phylogenies[ 0 ].setName( tree_name );
            }
            if ( !ForesterUtil.isEmpty( tree_id ) ) {
                final String[] s_ary = tree_id.split( ":" );
                phylogenies[ 0 ].setIdentifier( new Identifier( s_ary[ 1 ], s_ary[ 0 ] ) );
            }
            if ( !ForesterUtil.isEmpty( tree_desc ) ) {
                for( int i = 0; i < phylogenies.length; ++i ) {
                    phylogenies[ i ].setDescription( tree_desc );
                }
            }
        }
        try {
            if ( advanced_table ) {
                Map<String, Map<String, String>> table = null;
                try {
                    table = PhylogenyDecorator.parseMappingTable( mapping_infile );
                }
                catch ( final IOException e ) {
                    ForesterUtil.fatalError( decorator.PRG_NAME, "failed to read \"" + mapping_infile + "\" ["
                            + e.getMessage() + "]" );
                }
                PhylogenyDecorator.decorate( phylogenies,
                                             table,
                                             picky,
                                             numbers_of_chars_allowed_to_remove_if_not_found_in_map );
            }
            else {
                PhylogenyDecorator.decorate( phylogenies,
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
        catch ( final NullPointerException e ) {
            ForesterUtil.unexpectedFatalError( decorator.PRG_NAME, e );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( decorator.PRG_NAME, "failed to map [" + e + "]" );
        }
        try {
            final PhylogenyWriter w = new PhylogenyWriter();
            w.toPhyloXML( phylogenies, 0, phylogenies_outfile, ForesterUtil.getLineSeparator() );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( decorator.PRG_NAME, "failed to write output [" + e.getMessage() + "]" );
        }
        System.out.println();
        ForesterUtil.programMessage( PRG_NAME, "wrote: " + phylogenies_outfile );
        ForesterUtil.programMessage( PRG_NAME, "OK." );
    }
}
