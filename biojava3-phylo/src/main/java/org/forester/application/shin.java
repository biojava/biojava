// $Id: shin.java,v 1.2 2009/10/15 01:13:10 cmzmasek Exp $
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
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.sdi.Shin;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public class shin {

    final static private String HELP_OPTION_1   = "help";
    final static private String HELP_OPTION_2   = "h";
    final static private String DEFAULT_OUTFILE = "out";
    final static private String PRG_NAME        = "shin";
    final static private String PRG_VERSION     = "0.001 alpha";
    final static private String PRG_DATE        = "2009.10.14";

    public static void main( final String args[] ) {
        ForesterUtil.printProgramInformation( PRG_NAME, PRG_VERSION, PRG_DATE );
        System.out.println();
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
        else if ( ( args.length != 3 ) ) {
            System.out.println();
            System.out.println( "wrong number of arguments" );
            System.out.println();
            print_help();
            System.exit( -1 );
        }
        final List<String> allowed_options = new ArrayList<String>();
        final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
        if ( dissallowed_options.length() > 0 ) {
            ForesterUtil.fatalError( PRG_NAME, "unknown option(s): " + dissallowed_options );
        }
        File gene_trees_dir = null;
        File species_trees_file = null;
        //File out_file = null;
        File out_dir = null;
        Phylogeny[] species_trees = null;
        try {
            gene_trees_dir = cla.getFile( 0 );
            species_trees_file = cla.getFile( 1 );
            out_dir = cla.getFile( 2 );
        }
        catch ( final IllegalArgumentException e ) {
            ForesterUtil.fatalError( PRG_NAME, "error in command line: " + e.getMessage() );
        }
        if ( ForesterUtil.isReadableFile( species_trees_file ) != null ) {
            ForesterUtil.fatalError( PRG_NAME, ForesterUtil.isReadableFile( species_trees_file ) );
        }
        if ( !gene_trees_dir.isDirectory() || !gene_trees_dir.canRead() ) {
            ForesterUtil.fatalError( PRG_NAME, "cannot read gene trees from [" + gene_trees_dir + "]" );
        }
        // if ( ForesterUtil.isWritableFile( out_file ) != null ) {
        //     ForesterUtil.fatalError( PRG_NAME, ForesterUtil.isWritableFile( out_file ) );
        // }
        if ( !out_dir.exists() ) {
            boolean success = false;
            try {
                success = out_dir.mkdir();
            }
            catch ( final Exception e ) {
                ForesterUtil.fatalError( PRG_NAME, "failed to create [" + out_dir + "] [" + e.getMessage() + "]" );
            }
            if ( !success ) {
                ForesterUtil.fatalError( PRG_NAME, "failed to create [" + out_dir + "]" );
            }
        }
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            species_trees = factory.create( species_trees_file, new PhyloXmlParser() );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to read species trees from [" + species_trees_file + "] ["
                    + e.getMessage() + "]" );
        }
        if ( ( species_trees == null ) || ( species_trees.length < 1 ) ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to read species trees from [" + species_trees_file + "]" );
        }
        ForesterUtil.programMessage( PRG_NAME, "read in " + species_trees.length + " species trees from ["
                + species_trees_file + "]" );
        final FilenameFilter filter = new FilenameFilter() {

            public boolean accept( final File dir, final String name ) {
                return ( !name.startsWith( "." ) && !name.startsWith( "00_" ) && name.endsWith( ".xml" ) );
            }
        };
        final String[] gene_tree_names = gene_trees_dir.list( filter );
        Arrays.sort( gene_tree_names );
        final List<File> gene_tree_files = new ArrayList<File>();
        for( final String gene_tree_name : gene_tree_names ) {
            final File gene_tree_file = new File( gene_trees_dir + ForesterUtil.FILE_SEPARATOR + gene_tree_name );
            if ( !gene_tree_file.isDirectory() ) {
                gene_tree_files.add( gene_tree_file );
            }
        }
        ForesterUtil.programMessage( PRG_NAME, "going to analyze " + gene_tree_files.size() + " gene trees from ["
                + gene_trees_dir + "]" );
        final Shin shin = new Shin();
        try {
            shin.method1( gene_tree_files, species_trees, out_dir );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
            e.printStackTrace();
        }
        ForesterUtil.programMessage( PRG_NAME, "OK" );
        //        System.out.println();
        //        System.out.println( "Strip species tree: " + strip );
        //        SDI sdi = null;
        //        final long start_time = new Date().getTime();
        //        try {
        //            if ( use_sdise ) {
        //                System.out.println();
        //                System.out.println( "Using SDIse algorithm." );
        //                sdi = new SDIse( gene_tree, species_tree );
        //            }
        //            else {
        //                System.out.println();
        //                System.out.println( "Using GSDI algorithm." );
        //                System.out.println();
        //                System.out.println( "Use most parsimonous duplication model: " + most_parsimonous_duplication_model );
        //                sdi = new GSDI( gene_tree, species_tree, most_parsimonous_duplication_model );
        //            }
        //        }
        //        catch ( final Exception e ) {
        //            ForesterUtil.unexpectedFatalError( PRG_NAME, e );
        //        }
        //        System.out.println();
        //        System.out.println( "Running time (excluding I/O): " + ( new Date().getTime() - start_time ) + "ms" );
        //        try {
        //            final PhylogenyWriter writer = new PhylogenyWriter();
        //            writer.toPhyloXML( out_file, gene_tree, 1 );
        //        }
        //        catch ( final IOException e ) {
        //            ForesterUtil.fatalError( PRG_NAME, "Failed to write to \"" + out_file + "\" [" + e.getMessage() + "]" );
        //        }
        //        System.out.println();
        //        System.out.println( "Successfully wrote resulting gene tree to: " + out_file );
        //        System.out.println();
        //        System.out.println();
    }

    private static void print_help() {
        System.out.println( "Usage: " + PRG_NAME + " [-options] <gene trees dir> <species tree file name> <outdir>" );
        System.out.println();
        System.out.println( "Options:" );
        System.out.println();
    }
}
