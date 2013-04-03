// $Id: goac.java,v 1.13 2009/10/30 02:53:54 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2009 Christian M. Zmasek
// Copyright (C) 2009 Burnham Institute for Medical Research
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
import java.util.SortedMap;

import org.forester.go.GoId;
import org.forester.go.GoTerm;
import org.forester.go.GoUtils;
import org.forester.go.OBOparser;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.CommandLineArguments;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterUtil;
import org.forester.util.GeneralTable;

public class goac {

    private static final String ALL           = "{ALL}";
    final static private String HELP_OPTION_1 = "help";
    final static private String HELP_OPTION_2 = "h";
    final static private String PRG_NAME      = "goac";
    final static private String PRG_VERSION   = "1.03";
    final static private String PRG_DATE      = "2009.05.08";
    final static private String E_MAIL        = "czmasek@burnham.org";
    final static private String WWW           = "www.phylosoft.org/forester/goac";

    private static void addStats( final SortedMap<String, List<GoId>> data_to_be_analyzed,
                                  final GeneralTable<String, Double> table ) {
        for( final String go : table.getColumnIdentifiers() ) {
            final DescriptiveStatistics stats = new BasicDescriptiveStatistics();
            for( final String label : data_to_be_analyzed.keySet() ) {
                if ( !label.equals( ALL ) ) {
                    final Double value = table.getValue( go, label );
                    stats.addValue( value == null ? 0.0 : value );
                }
            }
            table.setValue( go, "{AVG}", stats.arithmeticMean() );
            table.setValue( go, "{SUM}", stats.getSum() );
            table.setValue( go, "{MED}", stats.median() );
            table.setValue( go, "{SD}", stats.sampleStandardDeviation() );
            table.setValue( go, "{MIN}", stats.getMin() );
            table.setValue( go, "{MAX}", stats.getMax() );
        }
    }

    public static void main( final String args[] ) {
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
        final List<String> allowed_options = new ArrayList<String>();
        if ( cla.getNumberOfNames() != 3 ) {
            printHelp();
            System.exit( -1 );
        }
        final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
        if ( dissallowed_options.length() > 0 ) {
            ForesterUtil.fatalError( PRG_NAME, "unknown option(s): " + dissallowed_options );
        }
        final File obofile = cla.getFile( 0 );
        final File query_superterms_file = cla.getFile( 1 );
        final File exp_file = cla.getFile( 2 );
        final OBOparser parser = new OBOparser( obofile, OBOparser.ReturnType.BASIC_GO_TERM );
        List<GoTerm> all_go_terms = null;
        try {
            all_go_terms = parser.parse();
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, e.toString() );
        }
        final Map<GoId, GoTerm> goid_to_term_map = GoUtils.createGoIdToGoTermMap( all_go_terms );
        final List<GoId> query_superterms_ids = new ArrayList<GoId>();
        SortedMap<String, List<GoId>> query_superterms_id_raw = null;
        try {
            query_superterms_id_raw = GoUtils.parseGoIds( query_superterms_file, "#", "" );
        }
        catch ( final IOException e ) {
            ForesterUtil.printErrorMessage( PRG_NAME, e.getMessage() );
            System.exit( -1 );
        }
        final List<GoId> queries = query_superterms_id_raw.get( "" );
        for( final GoId id : queries ) {
            if ( !goid_to_term_map.containsKey( id ) ) {
                ForesterUtil.printErrorMessage( PRG_NAME, "\"" + id + "\" not present in \"" + obofile + "\"" );
                System.exit( -1 );
            }
            query_superterms_ids.add( id );
        }
        SortedMap<String, List<GoId>> data_to_be_analyzed = null;
        try {
            data_to_be_analyzed = GoUtils.parseGoIds( exp_file, "#", ">" );
        }
        catch ( final IOException e ) {
            ForesterUtil.printErrorMessage( PRG_NAME, e.getMessage() );
            System.exit( -1 );
        }
        final List<GoId> all_ids = new ArrayList<GoId>();
        for( final String label : data_to_be_analyzed.keySet() ) {
            final List<GoId> experiment_set_ids = data_to_be_analyzed.get( label );
            for( final GoId go_id : experiment_set_ids ) {
                if ( !goid_to_term_map.containsKey( go_id ) ) {
                    ForesterUtil.printErrorMessage( PRG_NAME, "GO id [" + go_id + "] not found in GO id to term map" );
                    System.exit( -1 );
                }
                all_ids.add( go_id );
            }
        }
        if ( data_to_be_analyzed.size() > 1 ) {
            data_to_be_analyzed.put( ALL, all_ids );
        }
        final GeneralTable<String, Double> table_counts = new GeneralTable<String, Double>();
        final GeneralTable<String, Double> table_percentage = new GeneralTable<String, Double>();
        for( final String label : data_to_be_analyzed.keySet() ) {
            System.out.println();
            System.out.println( label + "\t\t\t\t" );
            final List<GoId> experiment_set_ids = data_to_be_analyzed.get( label );
            Map<GoId, Integer> counts_id = null;
            try {
                counts_id = GoUtils.countCategoriesId( query_superterms_ids, experiment_set_ids, goid_to_term_map );
            }
            catch ( final Exception e ) {
                ForesterUtil.printErrorMessage( PRG_NAME, e.getMessage() );
                System.exit( -1 );
            }
            int sum = 0;
            for( final GoId id : counts_id.keySet() ) {
                sum += counts_id.get( id );
            }
            if ( sum > 0 ) {
                table_counts.setValue( "{total}", label, ( double ) sum );
            }
            for( final GoId id : counts_id.keySet() ) {
                final int counts = counts_id.get( id );
                double percentage = 0.0;
                if ( sum > 0 ) {
                    percentage = 100.0 * counts / ( sum );
                }
                System.out.println( counts + "\t" + counts + "/" + sum + "\t" + percentage + "\t" + id + "\t"
                        + goid_to_term_map.get( id ).getName() );
                table_counts.setValue( goid_to_term_map.get( id ).getName(), label, ( double ) counts );
                table_percentage.setValue( goid_to_term_map.get( id ).getName(), label, percentage );
            }
        }
        addStats( data_to_be_analyzed, table_counts );
        addStats( data_to_be_analyzed, table_percentage );
        System.out.println();
        System.out.println();
        System.out.println();
        System.out.println( table_counts.toString( ForesterUtil.FORMATTER_3 ) );
        System.out.println();
        System.out.println();
        System.out.println();
        System.out.println( table_percentage.toString( ForesterUtil.FORMATTER_3 ) );
        System.out.println();
    }

    private static void printHelp() {
        ForesterUtil.printProgramInformation( PRG_NAME, PRG_VERSION, PRG_DATE, E_MAIL, WWW );
        System.out.println( "Usage:" );
        System.out.println();
        System.out
                .println( PRG_NAME
                        + " <file with all GO terms, in 'obo' format> <file with ancestral term ids> <file with go ids to be analyzed>" );
        System.out.println();
        System.out.println();
    }
}
