// $Id: nj.java,v 1.7 2008/03/11 00:29:26 cmzmasek Exp $
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
import java.util.Date;
import java.util.List;
import org.biojava3.phylo.CheckTreeAccuracy;

import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogenyinference.DistanceMatrix;
import org.forester.phylogenyinference.NeighborJoining;
import org.forester.phylogenyinference.SymmetricalDistanceMatrixParser;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public class nj {

    final static private String HELP_OPTION_1         = "help";
    final static private String HELP_OPTION_2         = "h";
    final static private String VERBOSE_OPTION        = "v";
    final static private String UPPER_TRIANGLE_OPTION = "u";
    final static private String PRG_NAME              = "nj";
    final static private String PRG_VERSION           = "0.0.1";
    final static private String PRG_DATE              = "2008.03.04";
    final static private String E_MAIL                = "czmasek@burnham.org";
    final static private String WWW                   = "www.phylosoft.org/forester/";

    public static void main(  String args[] ) {
        String[] args1 = {"/Users/Scooter/mutualinformation/project/nuclear_receptor/PF00104_small.fasta.phylip","/Users/Scooter/mutualinformation/project/nuclear_receptor/PF00104.small.forester.temp"};
        args = args1;
        ForesterUtil.printProgramInformation( PRG_NAME, PRG_VERSION, PRG_DATE, E_MAIL, WWW );
        final List<String> allowed_options = new ArrayList<String>();
        allowed_options.add( HELP_OPTION_1 );
        allowed_options.add( HELP_OPTION_2 );
        allowed_options.add( VERBOSE_OPTION );
        allowed_options.add( UPPER_TRIANGLE_OPTION );
        if ( ( args.length < 2 ) ) {
            printHelp();
            System.exit( -1 );
        }
        CommandLineArguments cla = null;
        try {
            cla = new CommandLineArguments( args );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        if ( cla.isOptionSet( HELP_OPTION_1 ) || cla.isOptionSet( HELP_OPTION_2 ) ) {
            printHelp();
            System.exit( 0 );
        }
        if ( cla.getNumberOfNames() != 2 ) {
            printHelp();
            System.exit( -1 );
        }
        boolean verbose = false;
        boolean upper_triangle = false;
        if ( cla.isOptionSet( VERBOSE_OPTION ) ) {
            verbose = true;
        }
        if ( cla.isOptionSet( UPPER_TRIANGLE_OPTION ) ) {
            upper_triangle = true;
        }
        verbose = true;
        final File infile = cla.getFile( 0 );
        final File outfile = cla.getFile( 1 );
        final String error1 = ForesterUtil.isReadableFile( infile );
        if ( !ForesterUtil.isEmpty( error1 ) ) {
            ForesterUtil.fatalError( PRG_NAME, "cannot read from infile [" + infile + "]: " + error1 );
        }
        if ( outfile.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "outfile [" + outfile + "] already exists" );
        }
        final String error2 = ForesterUtil.isWritableFile( outfile );
        if ( !ForesterUtil.isEmpty( error2 ) ) {
            ForesterUtil.fatalError( PRG_NAME, "cannot write to outfile [" + outfile + "]: " + error2 );
        }
        final SymmetricalDistanceMatrixParser parser = SymmetricalDistanceMatrixParser.createInstance();
        if ( upper_triangle ) {
            parser.setInputMatrixType( SymmetricalDistanceMatrixParser.InputMatrixType.UPPER_TRIANGLE );
        }
        else {
            parser.setInputMatrixType( SymmetricalDistanceMatrixParser.InputMatrixType.LOWER_TRIANGLE );
        }
        DistanceMatrix[] matrices = null;
        DistanceMatrix distanceMatrix = null;
        try {
            matrices = parser.parse( infile );            
            distanceMatrix = CheckTreeAccuracy.copyMatrix(matrices[0]);         
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to read from infile [" + infile + "]: " + e.getMessage() );
        }
        if ( verbose ) {
            System.out.println( PRG_NAME + " > read " + matrices.length + " pairwise distance matrice(s) of size "
                    + matrices[ 0 ].getSize() );
        }
        final List<Phylogeny> ps = new ArrayList<Phylogeny>();
        final NeighborJoining nj = NeighborJoining.createInstance();
        nj.setVerbose( verbose );
        final long start_time = new Date().getTime();
        for( final DistanceMatrix matrix : matrices ) {
            ps.add( nj.execute( matrix ) );
        }
        final long end_time = new Date().getTime();

        Phylogeny p = ps.get(0);
        CheckTreeAccuracy checkTreeAccuracy = new CheckTreeAccuracy();
        checkTreeAccuracy.process(p,distanceMatrix );

        final PhylogenyWriter w = new PhylogenyWriter();
        try {
            w.toPhyloXML( outfile, ps, 1, ForesterUtil.LINE_SEPARATOR );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to write to outfile [" + outfile + "]: " + e.getMessage() );
        }
        System.out.println();
        System.out.println( PRG_NAME + " > OK [" + ( end_time - start_time ) + "ms]" );
        System.out.println();
    }

    private static void printHelp() {
        System.out.println();
        System.out.println( "Usage:" );
        System.out.println();
        System.out.println( "% java -cp forester.jar org.forester.applications." + PRG_NAME
                + " [options] <pairwise distances infile> <out file>" );
        System.out.println();
        System.out.println( " Options: " );
        System.out.println( VERBOSE_OPTION + ": verbose on" );
        System.out.println( UPPER_TRIANGLE_OPTION + ": upper triangle option on (lower triangle is default)" );
        System.out.println();
    }
}
