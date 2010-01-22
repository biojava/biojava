// $Id: nhx_too.java,v 1.8 2009/11/20 22:22:09 cmzmasek Exp $
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
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public class nhx_too {

    final static private String PRG_NAME                 = "nhx_too";
    final static private String PRG_VERSION              = "0.1";
    final static private String PRG_DATE                 = "2008.03.04";
    final static private String INT_NODE_NAME_IS_SUPPORT = "is";

    public static void main( final String args[] ) {
        ForesterUtil.printProgramInformation( nhx_too.PRG_NAME, nhx_too.PRG_VERSION, nhx_too.PRG_DATE );
        if ( ( args.length < 3 ) || ( args.length > 3 ) ) {
            System.out.println();
            System.out.println( nhx_too.PRG_NAME + ": wrong number of arguments" );
            System.out.println();
            System.out.println( "Usage: \"" + nhx_too.PRG_NAME + " [options] <infile> <outfile>\n" );
            System.out.println( " Options: -" + nhx_too.INT_NODE_NAME_IS_SUPPORT
                    + ": internal node names are support values (i.e. MrBayes output)" );
            System.out.println();
            System.exit( -1 );
        }
        CommandLineArguments cla = null;
        try {
            cla = new CommandLineArguments( args );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        final List<String> allowed_options = new ArrayList<String>();
        allowed_options.add( nhx_too.INT_NODE_NAME_IS_SUPPORT );
        final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
        if ( dissallowed_options.length() > 0 ) {
            ForesterUtil.fatalError( nhx_too.PRG_NAME, "Unknown option(s): " + dissallowed_options );
        }
        final File phylogeny_infile = cla.getFile( 0 );
        final File phylogeny_outfile = cla.getFile( 1 );
        boolean int_node_name_is_support = false;
        if ( cla.isOptionSet( nhx_too.INT_NODE_NAME_IS_SUPPORT ) ) {
            int_node_name_is_support = true;
        }
        Phylogeny p = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ForesterUtil.createParserDependingOnFileType( phylogeny_infile, true );
            p = factory.create( phylogeny_infile, pp )[ 0 ];
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( nhx_too.PRG_NAME, "Could not read \"" + phylogeny_infile + "\" [" + e.getMessage()
                    + "]" );
        }
        if ( int_node_name_is_support ) {
            try {
                ForesterUtil.transferInternalNodeNamesToConfidence( p );
            }
            catch ( final Exception e ) {
                ForesterUtil.unexpectedFatalError( nhx_too.PRG_NAME,
                                                   "Failure during moving of internal names to support values ["
                                                           + e.getMessage() + "]" );
            }
        }
        try {
            final PhylogenyWriter w = new PhylogenyWriter();
            w.toNewHampshireX( p, phylogeny_outfile );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( nhx_too.PRG_NAME, "Failure to write output [" + e.getMessage() + "]" );
        }
        System.out.println();
        System.out.println( "Done [wrote \"" + phylogeny_outfile + "\"]." );
        System.out.println();
    }
}
