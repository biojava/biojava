// $Id: count_support.java,v 1.12 2009/11/20 22:22:10 cmzmasek Exp $
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
import java.util.Arrays;
import java.util.List;

import org.forester.development.SupportCount;
import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public class count_support {

    final static private String  PRG_NAME                = "count_support";
    final static private String  PRG_VERSION             = "1.0";
    final static private String  PRG_DATE                = "2008.03.04";
    private final static boolean WRITE_EVALUATORS_AS_NHX = false;

    public static void main( final String args[] ) {
        ForesterUtil
                .printProgramInformation( count_support.PRG_NAME, count_support.PRG_VERSION, count_support.PRG_DATE );
        if ( ( args.length < 3 ) || ( args.length > 7 ) ) {
            System.out.println();
            System.out.println( count_support.PRG_NAME + ": wrong number of arguments" );
            System.out.println();
            System.out
                    .println( "Usage: \"count_support [options] <file containing phylogeny to be evaluated> <file with phylogenies to be used for evaluation> <outfile> [outfile for evaluator phylogenies, "
                            + "always unstripped if -t=<d> option is used, otherwise strippedness is dependent on -s option]\"\n" );
            System.out
                    .println( " Options: -s strip external nodes from evaluator phylogenies not found in phylogeny to be evaluated" );
            System.out.println( "        : -t=<d> threshold for similarity (0.0 to 1.0)" );
            System.out.println( "        : -n no branch lengths in outfile for evaluator phylogenies" );
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
        allowed_options.add( "s" );
        allowed_options.add( "t" );
        allowed_options.add( "n" );
        final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
        if ( dissallowed_options.length() > 0 ) {
            ForesterUtil.fatalError( count_support.PRG_NAME, "Unknown option(s): " + dissallowed_options );
        }
        final File phylogeny_infile = cla.getFile( 0 );
        final File evaluators_infile = cla.getFile( 1 );
        final File phylogeny_outfile = cla.getFile( 2 );
        File evaluators_outfile = null;
        boolean branch_lengths_in_ev_out = true;
        if ( cla.isOptionSet( "n" ) ) {
            branch_lengths_in_ev_out = false;
        }
        if ( cla.getNumberOfNames() == 4 ) {
            evaluators_outfile = cla.getFile( 3 );
        }
        else {
            if ( !branch_lengths_in_ev_out ) {
                ForesterUtil.fatalError( count_support.PRG_NAME,
                                         "Cannot use -n option if no outfile for evaluators specified" );
            }
        }
        Phylogeny p = null;
        Phylogeny[] ev = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ForesterUtil.createParserDependingOnFileType( phylogeny_infile, true );
            p = factory.create( phylogeny_infile, pp )[ 0 ];
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( count_support.PRG_NAME, "Could not read \"" + phylogeny_infile + "\" ["
                    + e.getMessage() + "]" );
        }
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ForesterUtil.createParserDependingOnFileType( evaluators_infile, true );
            ev = factory.create( evaluators_infile, pp );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( count_support.PRG_NAME, "Could not read \"" + evaluators_infile + "\" ["
                    + e.getMessage() + "]" );
        }
        boolean strip = false;
        if ( cla.isOptionSet( "s" ) ) {
            strip = true;
        }
        double threshhold = -1.0;
        if ( cla.isOptionSet( "t" ) ) {
            try {
                threshhold = cla.getOptionValueAsDouble( "t" );
            }
            catch ( final Exception e ) {
                ForesterUtil.fatalError( count_support.PRG_NAME, "error in command line arguments: " + e.getMessage() );
            }
            if ( ( threshhold < 0 ) || ( threshhold > 1.0 ) ) {
                ForesterUtil.fatalError( count_support.PRG_NAME,
                                         "support threshold has to be between 0.0 and 1.0 (inclusive)" );
            }
        }
        List<Phylogeny> evaluator_phylogenies_above_threshold = null;
        try {
            if ( threshhold >= 0 ) {
                evaluator_phylogenies_above_threshold = SupportCount.count( p, ev, strip, threshhold, true );
                if ( evaluator_phylogenies_above_threshold.size() < 1 ) {
                    ForesterUtil.fatalError( "count_support", "appears like threshold for similarity is set too high" );
                }
            }
            else {
                SupportCount.count( p, ev, strip, true );
            }
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( count_support.PRG_NAME, "Failure during support counting: " + e.getMessage() );
        }
        if ( threshhold >= 0 ) {
            count_support.normalizeSupport( p, 100, evaluator_phylogenies_above_threshold.size() );
            System.out.println( evaluator_phylogenies_above_threshold.size() + " out of " + ev.length
                    + " evaluator phylogenies are above threshold of " + threshhold );
        }
        try {
            final PhylogenyWriter w = new PhylogenyWriter();
            w.toPhyloXML( phylogeny_outfile, p, 1 );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( count_support.PRG_NAME, "Failure to write output [" + e.getMessage() + "]" );
        }
        System.out.println();
        System.out.println( "Wrote phylogeny with support values to: " + phylogeny_outfile );
        if ( evaluators_outfile != null ) {
            try {
                final PhylogenyWriter w = new PhylogenyWriter();
                if ( evaluator_phylogenies_above_threshold != null ) {
                    System.out.println( "Writing " + evaluator_phylogenies_above_threshold.size()
                            + " evaluator phylogenies above threshold of " + threshhold + " to: " + evaluators_outfile );
                    if ( count_support.WRITE_EVALUATORS_AS_NHX ) {
                        w.toNewHampshireX( evaluator_phylogenies_above_threshold, evaluators_outfile, ";"
                                + ForesterUtil.getLineSeparator() );
                    }
                    else {
                        w.toNewHampshire( evaluator_phylogenies_above_threshold,
                                          true,
                                          branch_lengths_in_ev_out,
                                          evaluators_outfile,
                                          ";" + ForesterUtil.getLineSeparator() );
                    }
                }
                else {
                    System.out.println( "Writing " + ev.length + " evaluator phylogenies to :" + evaluators_outfile );
                    if ( count_support.WRITE_EVALUATORS_AS_NHX ) {
                        w.toNewHampshireX( Arrays.asList( ev ), evaluators_outfile, ";"
                                + ForesterUtil.getLineSeparator() );
                    }
                    else {
                        w.toNewHampshire( Arrays.asList( ev ), true, branch_lengths_in_ev_out, evaluators_outfile, ";"
                                + ForesterUtil.getLineSeparator() );
                    }
                }
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( count_support.PRG_NAME, "Failure to write output [" + e.getMessage() + "]" );
            }
        }
        System.out.println();
        System.out.println( "Done." );
        System.out.println();
    }

    private static void normalizeSupport( final Phylogeny p, final double normalized_max, final int number_phylos ) {
        double min = Double.MAX_VALUE;
        double max = -Double.MAX_VALUE;
        double sum = 0.0;
        int n = 0;
        for( final PhylogenyNodeIterator iter = p.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( !node.isRoot() && !node.isExternal() ) {
                final double b = PhylogenyMethods.getConfidenceValue( node );
                if ( b > max ) {
                    max = b;
                }
                if ( ( b >= 0 ) && ( b < min ) ) {
                    min = b;
                }
                sum += b;
                ++n;
            }
        }
        double av = sum / n;
        System.out.println( "Max support before normalization is    : " + max );
        System.out.println( "Min support before normalization is    : " + min );
        System.out.println( "Average support before normalization is: " + av + " (=" + sum + "/" + n + ")" );
        System.out.println( "Normalizing so that theoretical maximum support value is: " + normalized_max );
        System.out.println( "Number of phylogenies used in support analysis: " + number_phylos );
        final double f = normalized_max / number_phylos;
        min = Double.MAX_VALUE;
        max = -Double.MAX_VALUE;
        sum = 0.0;
        n = 0;
        for( final PhylogenyNodeIterator iter = p.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( node.isRoot() || node.isExternal() ) {
                PhylogenyMethods.setBootstrapConfidence( node, Confidence.CONFIDENCE_DEFAULT_VALUE );
            }
            else {
                double b = PhylogenyMethods.getConfidenceValue( node );
                b = f * b;
                PhylogenyMethods.setBootstrapConfidence( node, b );
                if ( b > max ) {
                    max = b;
                }
                if ( ( b >= 0 ) && ( b < min ) ) {
                    min = b;
                }
                sum += b;
                ++n;
            }
        }
        av = sum / n;
        System.out.println( "Max support after normalization is    : " + max );
        System.out.println( "Min support after normalization is    : " + min );
        System.out.println( "Average support after normalization is: " + av + " (=" + sum + "/" + n + ")" );
    }
}
