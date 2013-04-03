// $Id: strip.java,v 1.11 2009/11/20 22:22:09 cmzmasek Exp $
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

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.util.ForesterUtil;

public class strip {

    public static void main( final String args[] ) {
        if ( args.length < 4 ) {
            System.out.println( "\nstrip: Wrong number of arguments.\n" );
            System.out
                    .println( "Usage: \"strip <infile> <outfile> <options> [name1] [name2] ... OR [phylogenyfile]\"\n" );
            System.out.println( " Options: -k to keep listed nodes" );
            System.out.println( "          -r to remove listed nodes" );
            System.out.println( "          -kp to keep nodes found in [phylogenyfile]" );
            System.out.println( "          -rp to remove nodes found in [phylogenyfile]\n" );
            System.exit( -1 );
        }
        final File infile = new File( args[ 0 ] );
        final File outfile = new File( args[ 1 ] );
        final String options = args[ 2 ];
        Phylogeny p = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ForesterUtil.createParserDependingOnFileType( infile, true );
            p = factory.create( infile, pp )[ 0 ];
        }
        catch ( final Exception e ) {
            System.out.println( "\nCould not read \"" + infile + "\" [" + e.getMessage() + "]\n" );
            System.exit( -1 );
        }
        boolean keep = false;
        boolean from_p0 = false;
        if ( options.trim().toLowerCase().equals( "-k" ) ) {
            keep = true;
        }
        else if ( options.trim().toLowerCase().equals( "-kp" ) ) {
            keep = true;
            from_p0 = true;
        }
        else if ( options.trim().toLowerCase().equals( "-rp" ) ) {
            from_p0 = true;
        }
        else if ( !options.trim().toLowerCase().equals( "-r" ) ) {
            System.out.println( "\nUnknown option \"" + options + "\"\n" );
            System.exit( -1 );
        }
        String[] names = null;
        if ( from_p0 ) {
            names = strip.readInNamesFromPhylogeny( args[ 3 ] );
        }
        else {
            names = new String[ args.length - 3 ];
            for( int i = 0; i < args.length - 3; ++i ) {
                names[ i ] = args[ i + 3 ];
            }
        }
        if ( keep ) {
            PhylogenyMethods.deleteExternalNodesPositiveSelection( names, p );
        }
        else {
            PhylogenyMethods.deleteExternalNodesNegativeSelection( names, p );
        }
        try {
            final PhylogenyWriter w = new PhylogenyWriter();
            w.toPhyloXML( outfile, p, 1 );
        }
        catch ( final IOException e ) {
            System.out.println( "\nFailure to write output [" + e.getMessage() + "]\n" );
            System.exit( -1 );
        }
    }

    private static String[] readInNamesFromPhylogeny( final String file ) {
        Phylogeny p0 = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final File f = new File( file );
            final PhylogenyParser pp = ForesterUtil.createParserDependingOnFileType( f, true );
            p0 = factory.create( f, pp )[ 0 ];
        }
        catch ( final Exception e ) {
            System.out.println( "\nCould not read \"" + file + "\" [" + e.getMessage() + "]\n" );
            System.exit( -1 );
        }
        return p0.getAllExternalNodeNames();
    }
}
