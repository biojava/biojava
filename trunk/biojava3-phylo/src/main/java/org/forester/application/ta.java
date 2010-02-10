// $Id: ta.java,v 1.4 2009/11/20 22:22:09 cmzmasek Exp $
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
import java.util.Date;
import java.util.List;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.sdi.GSDI;
import org.forester.sdi.SDI;
import org.forester.sdi.SDIse;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public class ta {

    final static private String STRIP_OPTION             = "s";
    final static private String SDISE_OPTION             = "b";
    final static private String MOST_PARSIMONIOUS_OPTION = "m";
    final static private String HELP_OPTION_1            = "help";
    final static private String HELP_OPTION_2            = "h";
    final static private String DEFAULT_OUTFILE          = "sdi_out.xml";
    final static private String PRG_NAME                 = "sdi";
    final static private String PRG_VERSION              = "alpha 0.3";
    final static private String PRG_DATE                 = "2008.03.04";

    public static void main( final String args[] ) {
        ForesterUtil.printProgramInformation( PRG_NAME, PRG_VERSION, PRG_DATE );
        CommandLineArguments cla = null;
        try {
            cla = new CommandLineArguments( args );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        if ( cla.isOptionSet( HELP_OPTION_1 ) || cla.isOptionSet( HELP_OPTION_2 ) ) {
            System.out.println();
            print_help();
            System.exit( 0 );
        }
        else if ( ( args.length < 2 ) || ( cla.getNumberOfNames() < 2 ) || ( cla.getNumberOfNames() > 3 ) ) {
            System.out.println();
            System.out.println( "Wrong number of arguments." );
            System.out.println();
            print_help();
            System.exit( -1 );
        }
        final List<String> allowed_options = new ArrayList<String>();
        allowed_options.add( STRIP_OPTION );
        allowed_options.add( SDISE_OPTION );
        allowed_options.add( MOST_PARSIMONIOUS_OPTION );
        final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
        if ( dissallowed_options.length() > 0 ) {
            ForesterUtil.fatalError( PRG_NAME, "unknown option(s): " + dissallowed_options );
        }
        boolean use_sdise = false;
        boolean strip = false;
        boolean most_parsimonous_duplication_model = false;
        if ( cla.isOptionSet( STRIP_OPTION ) ) {
            strip = true;
        }
        if ( cla.isOptionSet( SDISE_OPTION ) ) {
            use_sdise = true;
        }
        if ( cla.isOptionSet( MOST_PARSIMONIOUS_OPTION ) ) {
            if ( use_sdise ) {
                ForesterUtil.fatalError( PRG_NAME, "Can only use most parsimonious duplication mode with GSDI" );
            }
            most_parsimonous_duplication_model = true;
        }
        Phylogeny species_tree = null;
        Phylogeny gene_tree = null;
        File gene_tree_file = null;
        File species_tree_file = null;
        File out_file = null;
        try {
            gene_tree_file = cla.getFile( 0 );
            species_tree_file = cla.getFile( 1 );
            if ( cla.getNumberOfNames() == 3 ) {
                out_file = cla.getFile( 2 );
            }
            else {
                out_file = new File( DEFAULT_OUTFILE );
            }
        }
        catch ( final IllegalArgumentException e ) {
            ForesterUtil.fatalError( PRG_NAME, "error in command line: " + e.getMessage() );
        }
        if ( ForesterUtil.isReadableFile( gene_tree_file ) != null ) {
            ForesterUtil.fatalError( PRG_NAME, ForesterUtil.isReadableFile( gene_tree_file ) );
        }
        if ( ForesterUtil.isReadableFile( species_tree_file ) != null ) {
            ForesterUtil.fatalError( PRG_NAME, ForesterUtil.isReadableFile( species_tree_file ) );
        }
        if ( ForesterUtil.isWritableFile( out_file ) != null ) {
            ForesterUtil.fatalError( PRG_NAME, ForesterUtil.isWritableFile( out_file ) );
        }
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ForesterUtil.createParserDependingOnFileType( species_tree_file, true );
            species_tree = factory.create( species_tree_file, pp )[ 0 ];
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "Failed to read species tree from \"" + gene_tree_file + "\" ["
                    + e.getMessage() + "]" );
        }
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ForesterUtil.createParserDependingOnFileType( gene_tree_file, true );
            gene_tree = factory.create( gene_tree_file, pp )[ 0 ];
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "Failed to read gene tree from \"" + gene_tree_file + "\" ["
                    + e.getMessage() + "]" );
        }
        gene_tree.setRooted( true );
        species_tree.setRooted( true );
        if ( !gene_tree.isCompletelyBinary() ) {
            ForesterUtil.fatalError( PRG_NAME, "gene tree is not completely binary." );
        }
        if ( use_sdise ) {
            if ( !species_tree.isCompletelyBinary() ) {
                ForesterUtil.fatalError( PRG_NAME, "species tree is not completely binary." );
            }
        }
        // For timing.
        // gene_tree = Helper.createBalancedTree( 10 );
        // species_tree = Helper.createBalancedTree( 13 );
        // species_tree = Helper.createUnbalancedTree( 1024 );
        // gene_tree = Helper.createUnbalancedTree( 8192 );
        // species_tree = gene_tree.copyTree();
        // gene_tree = species_tree.copyTree();
        // Helper.numberSpeciesInOrder( species_tree );
        // Helper.numberSpeciesInOrder( gene_tree );
        // Helper.randomizeSpecies( 1, 8192, gene_tree );
        // Helper.intervalNumberSpecies( gene_tree, 4096 );
        // Helper.numberSpeciesInDescOrder( gene_tree );
        System.out.println();
        System.out.println( "Strip species tree: " + strip );
        SDI sdi = null;
        final long start_time = new Date().getTime();
        try {
            if ( use_sdise ) {
                System.out.println();
                System.out.println( "Using SDIse algorithm." );
                sdi = new SDIse( gene_tree, species_tree );
            }
            else {
                System.out.println();
                System.out.println( "Using GSDI algorithm." );
                System.out.println();
                System.out.println( "Use most parsimonous duplication model: " + most_parsimonous_duplication_model );
                sdi = new GSDI( gene_tree, species_tree, most_parsimonous_duplication_model );
            }
        }
        catch ( final Exception e ) {
            ForesterUtil.unexpectedFatalError( PRG_NAME, e );
        }
        System.out.println();
        System.out.println( "Running time (excluding I/O): " + ( new Date().getTime() - start_time ) + "ms" );
        try {
            final PhylogenyWriter writer = new PhylogenyWriter();
            writer.toPhyloXML( out_file, gene_tree, 1 );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "Failed to write to \"" + out_file + "\" [" + e.getMessage() + "]" );
        }
        System.out.println();
        System.out.println( "Successfully wrote resulting gene tree to: " + out_file );
        System.out.println();
        //        if ( use_sdise ) {
        //            computeMappingCostL();
        //            System.out.println( "Mapping cost                    : " + computeMappingCostL() );
        //        }
        //        System.out.println( "Number of duplications          : " + getDuplicationsSum() );
        if ( !use_sdise && !most_parsimonous_duplication_model ) {
            System.out.println( "Number of potential duplications: "
                    + ( ( GSDI ) sdi ).getSpeciationOrDuplicationEventsSum() );
        }
        if ( !use_sdise ) {
            System.out.println( "Number speciations              : " + ( ( GSDI ) sdi ).getSpeciationsSum() );
        }
        System.out.println();
    } // main( final String args[] )

    private static void print_help() {
        System.out.println( "Usage: \"" + PRG_NAME
                + " [-options] <gene tree file name> <species tree file name>  [outfile name]\"" );
        System.out.println();
        System.out.println( "Options:" );
        System.out.println( " -" + STRIP_OPTION + ": to strip the species tree prior to duplication inference" );
        System.out.println( " -" + SDISE_OPTION
                + ": to use SDIse algorithm instead of GSDI algorithm (for binary trees only, faster)" );
        System.out.println( " -" + MOST_PARSIMONIOUS_OPTION + ": use most parimonious duplication model for GSDI: " );
        System.out.println( "     assign nodes as speciations which would otherwise be assiged" );
        System.out.println( "     as unknown because of polytomies in the species tree" );
        System.out.println();
        System.out.println( "Species tree file" );
        System.out.println( " In NHX format, with species names in species name fields unless -n option" );
        System.out.println( " is used." );
        System.out.println();
        System.out.println( "Gene tree file" );
        System.out.println( " In NHX format, with species names in species name fields and sequence names" );
        System.out.println( " in sequence name fields." );
        System.out.println();
        System.out.println( "!! WARNING: GSDI algorithm is under development, please use SDIse (-b) instead  !!" );
        System.out.println();
    } // print_help()
}
