// $Id: support_statistics.java,v 1.12 2009/11/20 22:22:09 cmzmasek Exp $
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
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.CommandLineArguments;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterUtil;

public final class support_statistics {

    final static private int    PLACES            = 2;
    final static private String HELP_OPTION       = "help";
    final static private String OUTPUTFILE_OPTION = "o";
    final static private String PRG_NAME          = "support_statistics";
    final static private String PRG_VERSION       = "1.0";
    final static private String PRG_DATE          = "2008.08.29";

    private static StringBuffer analyze( final File[] phylogenies_infiles, final Phylogeny[] phylogenies ) {
        final DescriptiveStatistics[] dss = new DescriptiveStatistics[ phylogenies.length ];
        for( int i = 0; i < phylogenies.length; i++ ) {
            dss[ i ] = new BasicDescriptiveStatistics();
            final Phylogeny p = phylogenies[ i ];
            for( final PhylogenyNodeIterator iter = p.iteratorPostorder(); iter.hasNext(); ) {
                final PhylogenyNode node = iter.next();
                if ( !node.isRoot() && !node.isExternal() ) {
                    double s = PhylogenyMethods.getConfidenceValue( node );
                    if ( s < 0.0 ) {
                        s = 0.0;
                    }
                    dss[ i ].addValue( s );
                }
            }
        }
        DescriptiveStatistics dss_comp = null;
        if ( dss.length > 2 ) {
            dss_comp = new BasicDescriptiveStatistics();
            for( final DescriptiveStatistics element : dss ) {
                dss_comp.addValue( element.arithmeticMean() );
            }
        }
        int max_length = 30;
        for( int i = 0; i < phylogenies.length; i++ ) {
            final int l = phylogenies_infiles[ i ].getName().length();
            if ( l > max_length ) {
                max_length = l;
            }
        }
        final StringBuffer sb = new StringBuffer();
        sb.append( "\t" + ForesterUtil.normalizeString( "name:", max_length, true, ' ' ) + "\t" );
        sb.append( "median:" + "\t" );
        sb.append( "mean:" + "\t" );
        sb.append( "sd:" + "\t" );
        sb.append( "min:" + "\t" );
        sb.append( "max:" + "\t" );
        sb.append( "n:" + "\t" );
        if ( dss_comp != null ) {
            sb.append( "\"z-score\":" );
        }
        sb.append( ForesterUtil.getLineSeparator() );
        for( int i = 0; i < phylogenies.length; i++ ) {
            sb.append( i + 1 + ":\t"
                    + ForesterUtil.normalizeString( phylogenies_infiles[ i ].getName(), max_length, true, ' ' ) + "\t" );
            sb.append( ForesterUtil.round( dss[ i ].median(), support_statistics.PLACES ) + "\t" );
            sb.append( ForesterUtil.round( dss[ i ].arithmeticMean(), support_statistics.PLACES ) + "\t" );
            try {
                sb.append( ForesterUtil.round( dss[ i ].sampleStandardDeviation(), support_statistics.PLACES ) + "\t" );
            }
            catch ( final ArithmeticException ex ) {
                sb.append( "n/a\t" );
            }
            sb.append( ForesterUtil.round( dss[ i ].getMin(), support_statistics.PLACES ) + "\t" );
            sb.append( ForesterUtil.round( dss[ i ].getMax(), support_statistics.PLACES ) + "\t" );
            sb.append( dss[ i ].getN() + "\t" );
            if ( dss_comp != null ) {
                final double z_score = dss_comp.sampleStandardUnit( dss[ i ].arithmeticMean() );
                sb.append( ForesterUtil.round( z_score, support_statistics.PLACES ) + "\t" );
            }
            sb.append( ForesterUtil.getLineSeparator() );
        }
        if ( dss_comp != null ) {
            sb.append( ForesterUtil.getLineSeparator() );
            sb.append( "\t" + ForesterUtil.normalizeString( "values for support means:", max_length, true, ' ' )
                    + "\t\t" );
            sb.append( ForesterUtil.round( dss_comp.arithmeticMean(), support_statistics.PLACES ) + "\t" );
            sb.append( ForesterUtil.round( dss_comp.sampleStandardDeviation(), support_statistics.PLACES ) + "\t" );
            sb.append( ForesterUtil.round( dss_comp.getMin(), support_statistics.PLACES ) + "\t" );
            sb.append( ForesterUtil.round( dss_comp.getMax(), support_statistics.PLACES ) + "\t" );
        }
        return sb;
    }

    public static void main( final String args[] ) {
        ForesterUtil.printProgramInformation( support_statistics.PRG_NAME,
                                              support_statistics.PRG_VERSION,
                                              support_statistics.PRG_DATE );
        if ( ( args.length < 1 ) ) {
            System.out.println();
            System.out.println( "wrong number of arguments" );
            System.out.println();
            support_statistics.printHelp();
            System.exit( -1 );
        }
        CommandLineArguments cla = null;
        try {
            cla = new CommandLineArguments( args );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        if ( cla.isOptionSet( support_statistics.HELP_OPTION ) ) {
            System.out.println();
            support_statistics.printHelp();
            System.exit( 0 );
        }
        final List<String> allowed_options = new ArrayList<String>();
        allowed_options.add( support_statistics.OUTPUTFILE_OPTION );
        final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
        if ( dissallowed_options.length() > 0 ) {
            ForesterUtil.fatalError( support_statistics.PRG_NAME, "unknown option(s): " + dissallowed_options );
        }
        final File[] phylogenies_infiles = new File[ cla.getNumberOfNames() ];
        for( int i = 0; i < phylogenies_infiles.length; ++i ) {
            phylogenies_infiles[ i ] = cla.getFile( i );
        }
        File outfile = null;
        if ( cla.isOptionSet( support_statistics.OUTPUTFILE_OPTION ) ) {
            try {
                outfile = new File( cla.getOptionValue( support_statistics.OUTPUTFILE_OPTION ) );
            }
            catch ( final IllegalArgumentException e ) {
                ForesterUtil.fatalError( support_statistics.PRG_NAME, "error in command line: " + e.getMessage() );
            }
            final String error = ForesterUtil.isWritableFile( outfile );
            if ( error != null ) {
                ForesterUtil.fatalError( support_statistics.PRG_NAME, error );
            }
        }
        final Phylogeny[] phylogenies = new Phylogeny[ phylogenies_infiles.length ];
        for( int i = 0; i < phylogenies_infiles.length; i++ ) {
            try {
                final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
                final PhylogenyParser pp = ForesterUtil
                        .createParserDependingOnFileType( phylogenies_infiles[ i ], true );
                phylogenies[ i ] = factory.create( phylogenies_infiles[ i ], pp )[ 0 ];
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( support_statistics.PRG_NAME, "could not read \"" + phylogenies_infiles[ i ]
                        + "\": " + e.getMessage() );
            }
        }
        final StringBuffer sb = support_statistics.analyze( phylogenies_infiles, phylogenies );
        System.out.println();
        System.out.println( sb );
        System.out.println();
        if ( outfile != null ) {
            try {
                final PrintWriter out = new PrintWriter( outfile );
                out.println( sb );
                out.flush();
                out.close();
                System.out.println( "wrote file: " + outfile );
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( support_statistics.PRG_NAME, "failed to write output: " + e.getMessage() );
            }
        }
        System.out.println( support_statistics.PRG_NAME + ": successfully completed" );
        System.out.println();
    }

    private static void printHelp() {
        System.out.println( "usage:" );
        System.out.println();
        System.out.println( support_statistics.PRG_NAME + " [-o=<outfile>] <phylogeny infile 1> "
                + "<phylogeny infile 2> <phylogeny infile 3> ..." );
        System.out.println();
        System.out.println( " options: " );
        System.out.println();
        System.out.println( " -o=<outfile> : write output to file" );
        System.out.println();
    }
}
