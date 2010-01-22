// $Id: phyloxml_converter.java,v 1.18 2009/11/20 22:22:10 cmzmasek Exp $
//
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

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.nexus.NexusPhylogeniesParser;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;
import org.forester.util.ForesterUtil.PhylogenyNodeField;

public class phyloxml_converter {

    final static private String  HELP_OPTION_1                    = "help";
    final static private String  HELP_OPTION_2                    = "h";
    final static private String  FIELD_OPTION                     = "f";
    final static private String  FIELD_CLADE_NAME                 = "nn";
    final static private String  FIELD_TAXONOMY_CODE              = "tc";
    final static private String  FIELD_TAXONOMY_SCI_NAME          = "sn";
    final static private String  FIELD_TAXONOMY_COMM_NAME         = "cn";
    final static private String  FIELD_SEQUENCE_GENE_NAME         = "gn";
    final static private String  FIELD_SEQUENCE_SYMBOL            = "sy";
    final static private String  FIELD_DUMMY                      = "dummy";
    final static private String  INTERNAL_NAMES_ARE_BOOT_SUPPPORT = "i";
    final static private String  MIDPOINT_REROOT                  = "m";
    final static private String  EXTRACT_TAXONOMY                 = "xt";
    final static private String  EXTRACT_TAXONOMY_PF              = "xp";
    final static private String  ORDER_SUBTREES                   = "o";
    final static private String  NO_TREE_LEVEL_INDENDATION        = "ni";
    final static private String  REPLACE_UNDER_SCORES             = "ru";
    final static private String  ALLOW_QUOTATION_MARKS_AND_SPACES = "qs";
    final static private String  PRG_NAME                         = "phyloxml_converter";
    final static private String  PRG_VERSION                      = "1.20";
    final static private String  PRG_DATE                         = "2009.10.13";
    final static private String  E_MAIL                           = "czmasek@burnham.org";
    final static private String  WWW                              = "www.phylosoft.org/forester/";
    final static private boolean SPECIAL                          = false;

    public static void main( final String args[] ) {
        ForesterUtil.printProgramInformation( PRG_NAME, PRG_VERSION, PRG_DATE, E_MAIL, WWW );
        CommandLineArguments cla = null;
        try {
            cla = new CommandLineArguments( args );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        if ( cla.isOptionSet( HELP_OPTION_1 ) || cla.isOptionSet( HELP_OPTION_2 ) || ( args.length == 0 ) ) {
            printHelp();
            System.exit( 0 );
        }
        if ( args.length < 3 ) {
            System.out.println();
            System.out.println( "[" + PRG_NAME + "] incorrect number of arguments" );
            System.out.println();
            printHelp();
            System.exit( -1 );
        }
        final List<String> allowed_options = new ArrayList<String>();
        allowed_options.add( NO_TREE_LEVEL_INDENDATION );
        allowed_options.add( FIELD_OPTION );
        allowed_options.add( MIDPOINT_REROOT );
        allowed_options.add( ORDER_SUBTREES );
        allowed_options.add( INTERNAL_NAMES_ARE_BOOT_SUPPPORT );
        allowed_options.add( REPLACE_UNDER_SCORES );
        allowed_options.add( EXTRACT_TAXONOMY );
        allowed_options.add( EXTRACT_TAXONOMY_PF );
        allowed_options.add( ALLOW_QUOTATION_MARKS_AND_SPACES );
        if ( cla.getNumberOfNames() != 2 ) {
            System.out.println();
            System.out.println( "[" + PRG_NAME + "] incorrect number of arguments" );
            System.out.println();
            printHelp();
            System.exit( -1 );
        }
        final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
        if ( dissallowed_options.length() > 0 ) {
            ForesterUtil.fatalError( PRG_NAME, "unknown option(s): " + dissallowed_options );
        }
        final List<String> mandatory_options = new ArrayList<String>();
        mandatory_options.add( FIELD_OPTION );
        final String missing_options = cla.validateMandatoryOptionsAsString( mandatory_options );
        if ( missing_options.length() > 0 ) {
            ForesterUtil.fatalError( PRG_NAME, "missing option(s): " + missing_options );
        }
        if ( !cla.isOptionValueSet( FIELD_OPTION ) ) {
            System.out.println();
            printHelp();
            System.exit( -1 );
        }
        final String field_option_value = cla.getOptionValue( FIELD_OPTION );
        PhylogenyNodeField field = null;
        if ( field_option_value.equals( FIELD_CLADE_NAME ) ) {
            field = PhylogenyNodeField.CLADE_NAME;
        }
        else if ( field_option_value.equals( FIELD_TAXONOMY_CODE ) ) {
            field = PhylogenyNodeField.TAXONOMY_CODE;
        }
        else if ( field_option_value.equals( FIELD_TAXONOMY_SCI_NAME ) ) {
            field = PhylogenyNodeField.TAXONOMY_SCIENTIFIC_NAME;
        }
        else if ( field_option_value.equals( FIELD_TAXONOMY_COMM_NAME ) ) {
            field = PhylogenyNodeField.TAXONOMY_COMMON_NAME;
        }
        else if ( field_option_value.equals( FIELD_SEQUENCE_GENE_NAME ) ) {
            field = PhylogenyNodeField.SEQUENCE_NAME;
        }
        else if ( field_option_value.equals( FIELD_SEQUENCE_SYMBOL ) ) {
            field = PhylogenyNodeField.SEQUENCE_SYMBOL;
        }
        else if ( field_option_value.equals( FIELD_DUMMY ) ) {
        }
        else {
            ForesterUtil.fatalError( PRG_NAME, "unknown value for -\"" + FIELD_OPTION + "\" option: \""
                    + field_option_value + "\"" );
        }
        boolean int_values_are_boots = false;
        if ( cla.isOptionSet( INTERNAL_NAMES_ARE_BOOT_SUPPPORT ) ) {
            int_values_are_boots = true;
        }
        boolean midpoint_reroot = false;
        if ( cla.isOptionSet( MIDPOINT_REROOT ) ) {
            midpoint_reroot = true;
        }
        boolean order_subtrees = false;
        if ( cla.isOptionSet( ORDER_SUBTREES ) ) {
            order_subtrees = true;
        }
        boolean replace_underscores = false;
        if ( cla.isOptionSet( REPLACE_UNDER_SCORES ) ) {
            replace_underscores = true;
        }
        boolean allow_quotes = false;
        if ( cla.isOptionSet( ALLOW_QUOTATION_MARKS_AND_SPACES ) ) {
            allow_quotes = true;
        }
        boolean no_indendation = false;
        if ( cla.isOptionSet( NO_TREE_LEVEL_INDENDATION ) ) {
            no_indendation = true;
        }
        boolean extr_taxonomy = false;
        if ( cla.isOptionSet( EXTRACT_TAXONOMY ) ) {
            extr_taxonomy = true;
        }
        boolean extr_taxonomy_pf_only = false;
        if ( cla.isOptionSet( EXTRACT_TAXONOMY_PF ) ) {
            extr_taxonomy_pf_only = true;
        }
        final File infile = cla.getFile( 0 );
        final File outfile = cla.getFile( 1 );
        if ( outfile.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + outfile + "] already exists" );
        }
        if ( !infile.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + infile + "] does not exist" );
        }
        Phylogeny[] phys = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser parser = ForesterUtil.createParserDependingOnFileType( infile, true );
            if ( parser instanceof NHXParser ) {
                if ( ( field != PhylogenyNodeField.TAXONOMY_CODE )
                        && ( field != PhylogenyNodeField.TAXONOMY_COMMON_NAME )
                        && ( field != PhylogenyNodeField.TAXONOMY_SCIENTIFIC_NAME ) ) {
                    if ( extr_taxonomy_pf_only ) {
                        ( ( NHXParser ) parser ).setTaxonomyExtraction( TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
                        replace_underscores = false;
                    }
                    else if ( extr_taxonomy ) {
                        ( ( NHXParser ) parser ).setTaxonomyExtraction( TAXONOMY_EXTRACTION.YES );
                        replace_underscores = false;
                    }
                }
                else {
                    ( ( NHXParser ) parser ).setTaxonomyExtraction( TAXONOMY_EXTRACTION.NO );
                }
                ( ( NHXParser ) parser ).setReplaceUnderscores( replace_underscores );
                ( ( NHXParser ) parser ).setIgnoreQuotes( !allow_quotes );
            }
            else if ( parser instanceof NexusPhylogeniesParser ) {
                ( ( NexusPhylogeniesParser ) parser ).setReplaceUnderscores( replace_underscores );
                ( ( NexusPhylogeniesParser ) parser ).setIgnoreQuotes( !allow_quotes );
            }
            phys = factory.create( infile, parser );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to read phylogeny from [" + infile + "]: " + e.getMessage() );
        }
        if ( SPECIAL ) {
            for( final Phylogeny phy : phys ) {
                performSpecialProcessing( phy );
            }
        }
        if ( int_values_are_boots ) {
            for( final Phylogeny phy : phys ) {
                ForesterUtil.transferInternalNamesToBootstrapSupport( phy );
            }
        }
        if ( field != null ) {
            for( final Phylogeny phy : phys ) {
                ForesterUtil.transferNodeNameToField( phy, field );
            }
        }
        if ( midpoint_reroot ) {
            try {
                for( final Phylogeny phy : phys ) {
                    PhylogenyMethods.midpointRoot( phy );
                }
            }
            catch ( final Exception e ) {
                System.out.println( "" );
                ForesterUtil.printWarningMessage( PRG_NAME, "midpoint rerooting failed: " + e.getLocalizedMessage() );
            }
        }
        if ( order_subtrees ) {
            for( final Phylogeny phy : phys ) {
                phy.orderAppearance( true );
            }
        }
        try {
            final PhylogenyWriter writer = new PhylogenyWriter();
            if ( no_indendation ) {
                writer.setIndentPhyloxml( false );
            }
            writer.toPhyloXML( phys, 0, outfile, ForesterUtil.LINE_SEPARATOR );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to write to [" + outfile + "]: " + e.getMessage() );
        }
        System.out.println( "[" + PRG_NAME + "] wrote: [" + outfile + "]" );
        System.out.println( "[" + PRG_NAME + "] OK" );
        System.out.println();
    }

    private static void performSpecialProcessing( final Phylogeny phy ) {
        // Can place some kind of custom processing here.
        //        final Set<PhylogenyNode> remove_us = new HashSet<PhylogenyNode>();
        //        int counter = 0;
        //        for( final PhylogenyNodeIterator it = phy.iteratorPostorder(); it.hasNext(); ) {
        //            final PhylogenyNode node = it.next();
        //            final String name = node.getNodeName().toLowerCase();
        //            if ( name.startsWith( "environmental_samples" ) || name.startsWith( "unclassified" )
        //                    || name.startsWith( "bacteria" ) || name.startsWith( "other" )
        //                    || name.startsWith( "viroids" ) || name.startsWith( "viruses" ) ) {
        //                remove_us.add( node );
        //                System.out.println( counter++ );
        //            }
        //        }
        //        phy.hashIDs();
        //        for( final PhylogenyNode node : remove_us ) {
        //            if ( phy.getNode( node.getNodeId() ) != null ) {
        //                phy.deleteSubtree( node );
        //                System.out.println( "deleted: " + node );
        //            }
        //        }
        //        phy.hashIDs();
        //
        //        for( final PhylogenyNodeIterator it = phy.iteratorPostorder(); it.hasNext(); ) {
        //            final PhylogenyNode node = it.next();
        //            node.getNodeData().setTaxonomy( null );
        //        }
        //        phy.reRoot( phy.getFirstExternalNode() );
        //        PhylogenyMethods.midpointRoot( phy );
        //        phy.orderAppearance( true );
        for( final PhylogenyNodeIterator it = phy.iteratorPostorder(); it.hasNext(); ) {
            final PhylogenyNode node = it.next();
            final String name = node.getNodeName();
            if ( !ForesterUtil.isEmpty( name ) ) {
                //                final Taxonomy taxo = new Taxonomy();
                //                if ( node.isExternal() ) {
                //                    taxo.setTaxonomyCode( name );
                //                    node.getNodeData().setTaxonomy( taxo );
                //                }
                //                else if ( name.indexOf( '_' ) == -1 || name.length() > 6 ) {
                //                    taxo.setScientificName( name );
                //                    node.getNodeData().setTaxonomy( taxo );
                //                }
                //                node.setName( "" );
                //                if ( name.indexOf( "BF" ) >= 0 ) {
                //                    taxo.setTaxonomyCode( "BACFR" );
                //                }
                //                else if ( name.indexOf( "BT" ) >= 0 ) {
                //                    taxo.setTaxonomyCode( "BACTN" );
                //                }
                //                else if ( name.indexOf( "MXAN" ) >= 0 ) {
                //                    taxo.setTaxonomyCode( "MYXXD" );
                //                }
                //                else if ( name.indexOf( "STIAU" ) >= 0 ) {
                //                    taxo.setTaxonomyCode( "STIAU" );
                //                }
                //                else if ( name.indexOf( "BOVA" ) >= 0 ) {
                //                    taxo.setTaxonomyCode( "BACOV" );
                //                }
                //                else if ( name.indexOf( "BUNI" ) >= 0 ) {
                //                    taxo.setTaxonomyCode( "BACUN" );
                //                }
                //                else if ( name.indexOf( "Pgin" ) >= 0 ) {
                //                    taxo.setTaxonomyCode( "PORGI" );
                //                }
                //                else if ( name.equals( "3CGH" ) || name.equals( "3CK7" ) ) {
                //                    taxo.setTaxonomyCode( "BACTN" );
                //                }
                // node.getNodeData().setTaxonomy( taxo );
            }
        }
    }

    private static void printHelp() {
        System.out.println( "Usage:" );
        System.out.println();
        System.out
                .println( PRG_NAME
                        + " -"
                        + FIELD_OPTION
                        + "=<field option> [options] <infile in New Hamphshire, NHX, Nexus, ToL XML, or phyloXML format> <outfile>" );
        System.out.println();
        System.out.println( " field options: " );
        System.out.println();
        System.out.println( "   " + FIELD_CLADE_NAME + ": transfer name to node/clade name" );
        System.out.println( "   " + FIELD_TAXONOMY_CODE + ": transfer name to taxonomy code" );
        System.out.println( "   " + FIELD_TAXONOMY_SCI_NAME + ": transfer name to taxonomy scientific name" );
        System.out.println( "   " + FIELD_TAXONOMY_COMM_NAME + ": transfer name to taxonomy common name" );
        System.out.println( "   " + FIELD_SEQUENCE_GENE_NAME + ": transfer name to sequence name" );
        System.out.println( "   " + FIELD_SEQUENCE_SYMBOL + ": transfer name to sequence symbol" );
        System.out.println();
        System.out.println( " options: " );
        System.out.println( " -" + INTERNAL_NAMES_ARE_BOOT_SUPPPORT
                + " : internal names in NH or NHX tree are bootstrap support values" );
        System.out.println( " -" + REPLACE_UNDER_SCORES + ": replace all underscores with spaces" );
        System.out.println( " -" + ALLOW_QUOTATION_MARKS_AND_SPACES
                + ": allow quotation marks and spaces (default: ignore)" );
        System.out.println( " -" + MIDPOINT_REROOT + " : midpoint reroot" );
        System.out.println( " -" + ORDER_SUBTREES + " : order subtrees" );
        System.out
                .println( " -"
                        + EXTRACT_TAXONOMY
                        + ": extract taxonomy to taxonomy code from \"seqname_TAXON\"-style names (cannot be used with the following field options: "
                        + FIELD_TAXONOMY_CODE + ", " + FIELD_TAXONOMY_COMM_NAME + ", " + FIELD_TAXONOMY_SCI_NAME + ")" );
        System.out
                .println( " -"
                        + EXTRACT_TAXONOMY_PF
                        + ": extract taxonomy to taxonomy code from Pfam (\"seqname_TAXON/x-y\") style names only (cannot be used with the following field options: "
                        + FIELD_TAXONOMY_CODE + ", " + FIELD_TAXONOMY_COMM_NAME + ", " + FIELD_TAXONOMY_SCI_NAME + ")" );
        System.out.println( " -" + NO_TREE_LEVEL_INDENDATION + ": no tree level indendation in phyloXML output" );
        System.out.println();
    }
}
