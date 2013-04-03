// $Id: MetaOntologizer.java,v 1.21 2009/10/26 23:29:40 cmzmasek Exp $
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
// WWW: www.phylosoft.org

package org.forester.go.etc;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Writer;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.go.GoId;
import org.forester.go.GoNameSpace;
import org.forester.go.GoTerm;
import org.forester.go.GoUtils;
import org.forester.go.OBOparser;
import org.forester.go.PfamToGoMapping;
import org.forester.surfacing.BasicSpecies;
import org.forester.surfacing.DomainId;
import org.forester.surfacing.Species;
import org.forester.surfacing.SurfacingConstants;
import org.forester.surfacing.SurfacingUtil;
import org.forester.util.ForesterUtil;

public class MetaOntologizer {

    private final static NumberFormat FORMATER                         = new DecimalFormat( "0.00E0" );
    private final static Color        MIN_COLOR                        = new Color( 0, 200, 50 );
    private final static Color        MAX_COLOR                        = new Color( 0, 0, 0 );
    final static private String       PRG_NAME                         = "meta_ontologizer";
    private static final boolean      VERBOSE                          = true;
    //table-a_41_dollo_all_gains_d-Topology-Elim-Bonferroni.txt:
    private final static Pattern      PATTERN_ONTOLOGIZER_TABLE_OUTPUT = Pattern
                                                                               .compile( ".*table-(.+?)_(dollo|fitch).*",
                                                                                         Pattern.CASE_INSENSITIVE );

    private static boolean hasResultsForSpecies( final Map<GoId, GoTerm> go_id_to_terms,
                                                 final SortedMap<String, SortedSet<OntologizerResult>> species_to_results_map,
                                                 final String species,
                                                 final GoNameSpace.GoNamespaceType namespace ) {
        for( final OntologizerResult ontologizer_result : species_to_results_map.get( species ) ) {
            if ( go_id_to_terms.get( ontologizer_result.getGoId() ).getGoNameSpace().getType() == namespace ) {
                return true;
            }
        }
        return false;
    }

    private static StringBuilder obtainDomainsForGoId( final List<PfamToGoMapping> pfam_to_go,
                                                       final SortedSet<DomainId> domains_per_species,
                                                       final Map<GoId, GoTerm> all_go_terms,
                                                       final GoId query_go_id,
                                                       final Set<DomainId> found_domain_ids ) {
        final StringBuilder sb = new StringBuilder();
        D: for( final DomainId domain_id : domains_per_species ) {
            for( final PfamToGoMapping ptg : pfam_to_go ) {
                if ( ptg.getKey().equals( domain_id ) ) {
                    final GoId go_id = ptg.getValue();
                    final Set<GoId> super_ids = new HashSet<GoId>();
                    for( final GoTerm term : GoUtils.getAllSuperGoTerms( go_id, all_go_terms ) ) {
                        super_ids.add( term.getGoId() );
                    }
                    super_ids.add( go_id );
                    if ( super_ids.contains( query_go_id ) ) {
                        sb.append( "[<a href=\"" + SurfacingConstants.PFAM_FAMILY_ID_LINK + domain_id + "\">"
                                + domain_id + "</a>] " );
                        found_domain_ids.add( domain_id );
                        continue D;
                    }
                }
            }
        }
        return sb;
    }

    private static String obtainSpecies( final File ontologizer_outfile ) {
        final Matcher matcher = PATTERN_ONTOLOGIZER_TABLE_OUTPUT.matcher( ontologizer_outfile.getName() );
        String species = null;
        if ( matcher.matches() ) {
            species = matcher.group( 1 );
            if ( VERBOSE ) {
                ForesterUtil
                        .programMessage( PRG_NAME, "species for [" + ontologizer_outfile + "] is [" + species + "]" );
            }
        }
        else {
            throw new IllegalStateException( "pattern [" + PATTERN_ONTOLOGIZER_TABLE_OUTPUT + "] did not match ["
                    + ontologizer_outfile.getName() + "]" );
        }
        return species;
    }

    private static SortedMap<Species, SortedSet<DomainId>> parseDomainGainLossFile( final File input )
            throws IOException {
        final String error = ForesterUtil.isReadableFile( input );
        if ( !ForesterUtil.isEmpty( error ) ) {
            throw new IOException( error );
        }
        final SortedMap<Species, SortedSet<DomainId>> speciesto_to_domain_id = new TreeMap<Species, SortedSet<DomainId>>();
        final BufferedReader br = new BufferedReader( new FileReader( input ) );
        String line;
        int line_number = 0;
        Species current_species = null;
        try {
            while ( ( line = br.readLine() ) != null ) {
                line_number++;
                line = line.trim();
                if ( ( ForesterUtil.isEmpty( line ) ) || ( line.startsWith( "##" ) ) ) {
                    // Ignore.
                }
                else if ( line.startsWith( "#" ) ) {
                    current_species = new BasicSpecies( line.substring( 1 ) );
                    speciesto_to_domain_id.put( current_species, new TreeSet<DomainId>() );
                }
                else {
                    if ( current_species == null ) {
                        throw new IOException( "parsing problem [at line " + line_number + "] in [" + input + "]" );
                    }
                    speciesto_to_domain_id.get( current_species ).add( new DomainId( line ) );
                }
            }
        }
        catch ( final Exception e ) {
            throw new IOException( "parsing problem [at line " + line_number + "] in [" + input + "]: "
                    + e.getMessage() );
        }
        return speciesto_to_domain_id;
    }

    private static void processOneSpecies( final Map<GoId, GoTerm> go_id_to_terms,
                                           final Writer b_html_writer,
                                           final Writer b_tab_writer,
                                           final Writer c_html_writer,
                                           final Writer c_tab_writer,
                                           final Writer m_html_writer,
                                           final Writer m_tab_writer,
                                           final SortedMap<String, SortedSet<OntologizerResult>> species_to_results_map,
                                           final String species,
                                           final double p_adjusted_upper_limit,
                                           final SortedSet<DomainId> domains_per_species,
                                           final List<PfamToGoMapping> pfam_to_go,
                                           final Set<DomainId> domain_ids_with_go_annot ) throws IOException {
        final SortedSet<OntologizerResult> ontologizer_results = species_to_results_map.get( species );
        for( final OntologizerResult ontologizer_result : ontologizer_results ) {
            final GoTerm go_term = go_id_to_terms.get( ontologizer_result.getGoId() );
            Writer current_html_writer = b_html_writer;
            Writer current_tab_writer = b_tab_writer;
            switch ( go_term.getGoNameSpace().getType() ) {
                case CELLULAR_COMPONENT:
                    current_html_writer = c_html_writer;
                    current_tab_writer = c_tab_writer;
                    break;
                case MOLECULAR_FUNCTION:
                    current_html_writer = m_html_writer;
                    current_tab_writer = m_tab_writer;
                    break;
            }
            writeValuesToTabWriter( species, ontologizer_result, go_term, current_tab_writer );
            writeValuesToHtmlWriter( ontologizer_result,
                                     go_term,
                                     current_html_writer,
                                     p_adjusted_upper_limit,
                                     species,
                                     go_id_to_terms,
                                     domains_per_species,
                                     pfam_to_go,
                                     domain_ids_with_go_annot );
        }
    }

    public static void reformat( final File ontologizer_outdir,
                                 final String result_file_prefix,
                                 final File domain_gain_loss_file,
                                 final String outfile_base,
                                 final File obo_file,
                                 final double p_adjusted_upper_limit,
                                 final String comment,
                                 final List<PfamToGoMapping> pfam_to_go ) throws IOException {
        if ( !ontologizer_outdir.exists() ) {
            throw new IllegalArgumentException( "[" + ontologizer_outdir + "] does not exist" );
        }
        if ( !ontologizer_outdir.isDirectory() ) {
            throw new IllegalArgumentException( "[" + ontologizer_outdir + "] is not a directory" );
        }
        if ( !obo_file.exists() ) {
            throw new IllegalArgumentException( "[" + obo_file + "] does not exist" );
        }
        if ( ( p_adjusted_upper_limit < 0.0 ) || ( p_adjusted_upper_limit > 1.0 ) ) {
            throw new IllegalArgumentException( "adjusted P values limit [" + p_adjusted_upper_limit
                    + "] is out of range" );
        }
        SortedMap<Species, SortedSet<DomainId>> speciesto_to_domain_id = null;
        if ( domain_gain_loss_file != null ) {
            if ( !domain_gain_loss_file.exists() ) {
                throw new IllegalArgumentException( "[" + domain_gain_loss_file + "] does not exist" );
            }
            speciesto_to_domain_id = parseDomainGainLossFile( domain_gain_loss_file );
            if ( VERBOSE ) {
                ForesterUtil.programMessage( PRG_NAME, "parsed gain/loss domains for " + speciesto_to_domain_id.size()
                        + " species from [" + domain_gain_loss_file + "]" );
            }
        }
        final String[] children = ontologizer_outdir.list();
        final List<File> ontologizer_outfiles = new ArrayList<File>();
        if ( children == null ) {
            throw new IllegalArgumentException( "problem with [" + ontologizer_outdir + "]" );
        }
        else {
            for( final String filename : children ) {
                if ( filename.startsWith( result_file_prefix ) ) {
                    ontologizer_outfiles.add( new File( filename ) );
                }
            }
        }
        if ( VERBOSE ) {
            ForesterUtil.programMessage( PRG_NAME, "need to analyze " + ontologizer_outfiles.size()
                    + " Ontologizer outfiles from [" + ontologizer_outdir + "]" );
        }
        final OBOparser parser = new OBOparser( obo_file, OBOparser.ReturnType.BASIC_GO_TERM );
        final List<GoTerm> go_terms = parser.parse();
        if ( VERBOSE ) {
            ForesterUtil.programMessage( PRG_NAME, "parsed " + go_terms.size() + " GO terms from [" + obo_file + "]" );
        }
        final Map<GoId, GoTerm> go_id_to_terms = GoUtils.createGoIdToGoTermMap( go_terms );
        if ( go_id_to_terms.size() != go_terms.size() ) {
            throw new IllegalArgumentException( "GO terms with non-unique ids found" );
        }
        final String b_file_html = outfile_base + "_B.html";
        final String b_file_txt = outfile_base + "_B.txt";
        final String m_file_html = outfile_base + "_C.html";
        final String m_file_txt = outfile_base + "_C.txt";
        final String c_file_html = outfile_base + "_M.html";
        final String c_file_txt = outfile_base + "_M.txt";
        final Writer b_html_writer = ForesterUtil.createBufferedWriter( b_file_html );
        final Writer b_tab_writer = ForesterUtil.createBufferedWriter( b_file_txt );
        final Writer c_html_writer = ForesterUtil.createBufferedWriter( m_file_html );
        final Writer c_tab_writer = ForesterUtil.createBufferedWriter( m_file_txt );
        final Writer m_html_writer = ForesterUtil.createBufferedWriter( c_file_html );
        final Writer m_tab_writer = ForesterUtil.createBufferedWriter( c_file_txt );
        final SortedMap<String, SortedSet<OntologizerResult>> species_to_results_map = new TreeMap<String, SortedSet<OntologizerResult>>();
        for( final File ontologizer_outfile : ontologizer_outfiles ) {
            final String species = obtainSpecies( ontologizer_outfile );
            final List<OntologizerResult> ontologizer_results = OntologizerResult.parse( new File( ontologizer_outdir
                    + ForesterUtil.FILE_SEPARATOR + ontologizer_outfile ) );
            final SortedSet<OntologizerResult> filtered_ontologizer_results = new TreeSet<OntologizerResult>();
            for( final OntologizerResult ontologizer_result : ontologizer_results ) {
                if ( ontologizer_result.getPAdjusted() <= p_adjusted_upper_limit ) {
                    filtered_ontologizer_results.add( ontologizer_result );
                }
            }
            species_to_results_map.put( species, filtered_ontologizer_results );
        }
        writeLabelsToTabWriter( b_tab_writer );
        writeLabelsToTabWriter( c_tab_writer );
        writeLabelsToTabWriter( m_tab_writer );
        String domain_gain_loss_file_full_path_str = null;
        if ( domain_gain_loss_file != null ) {
            domain_gain_loss_file_full_path_str = domain_gain_loss_file.getAbsolutePath();
        }
        writeHtmlHeader( b_html_writer,
                         GoNameSpace.GoNamespaceType.BIOLOGICAL_PROCESS.toString() + " | Pmax = "
                                 + p_adjusted_upper_limit + " | " + comment,
                         ontologizer_outdir.getAbsolutePath(),
                         domain_gain_loss_file_full_path_str );
        writeHtmlHeader( c_html_writer,
                         GoNameSpace.GoNamespaceType.CELLULAR_COMPONENT.toString() + " | Pmax = "
                                 + p_adjusted_upper_limit + " | " + comment,
                         ontologizer_outdir.getAbsolutePath(),
                         domain_gain_loss_file_full_path_str );
        writeHtmlHeader( m_html_writer,
                         GoNameSpace.GoNamespaceType.MOLECULAR_FUNCTION.toString() + " | Pmax = "
                                 + p_adjusted_upper_limit + " | " + comment,
                         ontologizer_outdir.getAbsolutePath(),
                         domain_gain_loss_file_full_path_str );
        for( final String species : species_to_results_map.keySet() ) {
            if ( hasResultsForSpecies( go_id_to_terms,
                                       species_to_results_map,
                                       species,
                                       GoNameSpace.GoNamespaceType.BIOLOGICAL_PROCESS ) ) {
                writeHtmlSpecies( b_html_writer, species );
            }
            if ( hasResultsForSpecies( go_id_to_terms,
                                       species_to_results_map,
                                       species,
                                       GoNameSpace.GoNamespaceType.CELLULAR_COMPONENT ) ) {
                writeHtmlSpecies( c_html_writer, species );
            }
            if ( hasResultsForSpecies( go_id_to_terms,
                                       species_to_results_map,
                                       species,
                                       GoNameSpace.GoNamespaceType.MOLECULAR_FUNCTION ) ) {
                writeHtmlSpecies( m_html_writer, species );
            }
            SortedSet<DomainId> domains_per_species = null;
            if ( ( speciesto_to_domain_id != null ) && ( speciesto_to_domain_id.size() > 0 ) ) {
                domains_per_species = speciesto_to_domain_id.get( new BasicSpecies( species ) );
            }
            final Set<DomainId> domain_ids_with_go_annot = new HashSet<DomainId>();
            processOneSpecies( go_id_to_terms,
                               b_html_writer,
                               b_tab_writer,
                               c_html_writer,
                               c_tab_writer,
                               m_html_writer,
                               m_tab_writer,
                               species_to_results_map,
                               species,
                               p_adjusted_upper_limit,
                               domains_per_species,
                               pfam_to_go,
                               domain_ids_with_go_annot );
            if ( ( speciesto_to_domain_id != null ) && ( speciesto_to_domain_id.size() > 0 ) ) {
                if ( hasResultsForSpecies( go_id_to_terms,
                                           species_to_results_map,
                                           species,
                                           GoNameSpace.GoNamespaceType.BIOLOGICAL_PROCESS ) ) {
                    writeHtmlDomains( b_html_writer, domains_per_species, domain_ids_with_go_annot );
                }
                if ( hasResultsForSpecies( go_id_to_terms,
                                           species_to_results_map,
                                           species,
                                           GoNameSpace.GoNamespaceType.CELLULAR_COMPONENT ) ) {
                    writeHtmlDomains( c_html_writer, domains_per_species, domain_ids_with_go_annot );
                }
                if ( hasResultsForSpecies( go_id_to_terms,
                                           species_to_results_map,
                                           species,
                                           GoNameSpace.GoNamespaceType.MOLECULAR_FUNCTION ) ) {
                    writeHtmlDomains( m_html_writer, domains_per_species, domain_ids_with_go_annot );
                }
            }
        }
        writeHtmlEnd( b_html_writer );
        writeHtmlEnd( c_html_writer );
        writeHtmlEnd( m_html_writer );
        b_html_writer.close();
        b_tab_writer.close();
        c_html_writer.close();
        c_tab_writer.close();
        m_html_writer.close();
        m_tab_writer.close();
        if ( VERBOSE ) {
            ForesterUtil.programMessage( PRG_NAME, "successfully wrote biological process summary to [" + b_file_html
                    + "]" );
            ForesterUtil.programMessage( PRG_NAME, "successfully wrote biological process summary to [" + b_file_txt
                    + "]" );
            ForesterUtil.programMessage( PRG_NAME, "successfully wrote molecular function summary to [" + m_file_html
                    + "]" );
            ForesterUtil.programMessage( PRG_NAME, "successfully wrote molecular function summary to [" + m_file_txt
                    + "]" );
            ForesterUtil.programMessage( PRG_NAME, "successfully wrote cellular component summary to [" + c_file_html
                    + "]" );
            ForesterUtil.programMessage( PRG_NAME, "successfully wrote cellular component summary to [" + c_file_txt
                    + "]" );
        }
    }

    private static void writeHtmlDomains( final Writer writer,
                                          final SortedSet<DomainId> domains,
                                          final Set<DomainId> domain_ids_with_go_annot ) throws IOException {
        writer.write( "<tr>" );
        writer.write( "<td colspan=\"10\">" );
        for( final DomainId domain : domains ) {
            if ( !domain_ids_with_go_annot.contains( domain ) ) {
                writer.write( "[<a class=\"new_type\" href=\"" + SurfacingConstants.PFAM_FAMILY_ID_LINK + domain
                        + "\">" + domain + "</a>] " );
            }
        }
        writer.write( "</td>" );
        writer.write( "</tr>" );
        writer.write( ForesterUtil.LINE_SEPARATOR );
    }

    private static void writeHtmlEnd( final Writer writer ) throws IOException {
        writer.write( "</table>" );
        writer.write( "</body>" );
        writer.write( "</html>" );
    }

    private static void writeHtmlHeader( final Writer w,
                                         final String desc,
                                         final String ontologizer_outdir,
                                         final String domain_gain_loss_file ) throws IOException {
        w.write( "<head>" );
        w.write( "<title>" );
        w.write( desc );
        w.write( "</title>" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( "<style>" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( "a:visited { color : #F87217; text-decoration : none; }" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( "a:link { color : #F87217; text-decoration : none; }" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( "a:hover { color : #FFFFFF; background-color : #00FF00; text-decoration : none; }" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( "a:hover { color : #FFFFFF; background-color : #00FF00; text-decoration : none; }" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( "a.new_type:visited { font-size: 7pt; color : #808080; text-decoration : none; }" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( "a.new_type:link { font-size: 7pt; color : #505050; text-decoration : none; }" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w
                .write( "a.new_type:hover { font-size: 7pt; color : #000000; background-color : #FFFF00; text-decoration : none; }" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w
                .write( "a.new_type:hover { font-size: 7pt; color : #000000; background-color : #FFFF00; text-decoration : none; }" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( "td { text-align: left; vertical-align: top; font-family: Verdana, Arial, Helvetica; font-size: 8pt}" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w
                .write( "th { text-align: left; vertical-align: top; font-family: Verdana, Arial, Helvetica; font-size: 10pt; font-weight: bold }" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( "h1 { color : #000000; font-family: Verdana, Arial, Helvetica; font-size: 18pt; font-weight: bold }" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( "h2 { color : #000000; font-family: Verdana, Arial, Helvetica; font-size: 16pt; font-weight: bold }" );
        w
                .write( "h3 { margin-top: 12px;  margin-bottom: 0px; color : #000000; font-family: Verdana, Arial, Helvetica; font-size: 12pt; font-weight: bold }" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( "</style>" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( "</head>" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( "<body>" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( "<h2>" );
        w.write( "meta ontologizer" );
        w.write( "</h2>" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( "<h2>" );
        w.write( desc );
        w.write( "</h2>" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( "<table>" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( "<tr><th>" );
        w.write( "ontolgizer output directory analysed:" );
        w.write( "</th><td>" );
        w.write( ontologizer_outdir );
        w.write( "</td></tr>" );
        if ( !ForesterUtil.isEmpty( domain_gain_loss_file ) ) {
            w.write( ForesterUtil.LINE_SEPARATOR );
            w.write( "<tr><th>" );
            w.write( "domain gain or loss file:" );
            w.write( "</th><td>" );
            w.write( domain_gain_loss_file );
            w.write( "</td></tr>" );
        }
        w.write( "</table>" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( "<table>" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( "<tr>" );
        w.write( "<th>" );
        w.write( "GO term name" );
        w.write( "</th><th>" );
        w.write( "GO id" );
        w.write( "</th><th>" );
        w.write( "P adjusted" );
        w.write( "</th><th>" );
        w.write( "P" );
        w.write( "</th><th>" );
        w.write( "Pop total" );
        w.write( "</th><th>" );
        w.write( "Pop term" );
        w.write( "</th><th>" );
        w.write( "Study total" );
        w.write( "</th><th>" );
        w.write( "Study term" );
        w.write( "</th><th>" );
        w.write( "Domains" );
        w.write( "</th><th>" );
        w.write( "trivial?" );
        w.write( "</th>" );
        w.write( "</tr>" );
        w.write( ForesterUtil.LINE_SEPARATOR );
    }

    private static void writeHtmlSpecies( final Writer writer, final String species ) throws IOException {
        writer.write( "<tr>" );
        writer.write( "<td><h3>" );
        writer.write( species );
        SurfacingUtil.writeTaxonomyLinks( writer, species );
        writer.write( "</h3></td>" );
        writer.write( "</tr>" );
        writer.write( ForesterUtil.LINE_SEPARATOR );
    }

    private static void writeLabelsToTabWriter( final Writer writer ) throws IOException {
        writer.write( "#species" );
        writer.write( "\t" );
        writer.write( "GO name" );
        writer.write( "\t" );
        writer.write( "GO id" );
        writer.write( "\t" );
        writer.write( "P adjusted" );
        writer.write( "\t" );
        writer.write( "P" );
        writer.write( "\t" );
        writer.write( "Pop total" );
        writer.write( "\t" );
        writer.write( "Pop term" );
        writer.write( "\t" );
        writer.write( "Study total" );
        writer.write( "\t" );
        writer.write( "Study term" );
        writer.write( "\t" );
        writer.write( "is trivial" );
        writer.write( ForesterUtil.LINE_SEPARATOR );
    }

    private static void writeValuesToHtmlWriter( final OntologizerResult ontologizer_result,
                                                 final GoTerm go_term,
                                                 final Writer writer,
                                                 final double p_adjusted_upper_limit,
                                                 final String species,
                                                 final Map<GoId, GoTerm> go_id_to_terms,
                                                 final SortedSet<DomainId> domains_per_species,
                                                 final List<PfamToGoMapping> pfam_to_go,
                                                 final Set<DomainId> domain_ids_with_go_annot ) throws IOException {
        final Color p_adj_color = ForesterUtil.calcColor( ontologizer_result.getPAdjusted(),
                                                          0,
                                                          p_adjusted_upper_limit,
                                                          MIN_COLOR,
                                                          MAX_COLOR );
        final Color p_color = ForesterUtil.calcColor( ontologizer_result.getP(),
                                                      0,
                                                      p_adjusted_upper_limit,
                                                      MIN_COLOR,
                                                      MAX_COLOR );
        writer.write( "<tr>" );
        writer.write( "<td>" );
        writer.write( "<font color=\"#" + ForesterUtil.colorToHex( p_adj_color ) + "\">" );
        writer.write( go_term.getName() );
        writer.write( "</font>" );
        writer.write( "</td><td>" );
        writer.write( "<a href=\"" + SurfacingConstants.GO_LINK + ontologizer_result.getGoId().getId()
                + "\" target=\"amigo_window\">" + ontologizer_result.getGoId().getId() + "</a>" );
        writer.write( "</td><td>" );
        writer.write( "<font color=\"#" + ForesterUtil.colorToHex( p_adj_color ) + "\">" );
        writer.write( FORMATER.format( ontologizer_result.getPAdjusted() ) );
        writer.write( "</font>" );
        writer.write( "</td><td>" );
        writer.write( "<font color=\"#" + ForesterUtil.colorToHex( p_color ) + "\">" );
        writer.write( FORMATER.format( ontologizer_result.getP() ) );
        writer.write( "</font>" );
        writer.write( "</td><td>" );
        writer.write( String.valueOf( ontologizer_result.getPopTotal() ) );
        writer.write( "</td><td>" );
        writer.write( String.valueOf( ontologizer_result.getPopTerm() ) );
        writer.write( "</td><td>" );
        writer.write( String.valueOf( ontologizer_result.getStudyTotal() ) );
        writer.write( "</td><td>" );
        writer.write( String.valueOf( ontologizer_result.getStudyTerm() ) );
        writer.write( "</td><td>" );
        if ( domains_per_species != null ) {
            final StringBuilder sb = obtainDomainsForGoId( pfam_to_go, domains_per_species, go_id_to_terms, go_term
                    .getGoId(), domain_ids_with_go_annot );
            writer.write( sb.toString() );
        }
        else {
            writer.write( " " );
        }
        writer.write( "</td><td>" );
        if ( ontologizer_result.isTrivial() ) {
            writer.write( "trivial" );
        }
        else {
            writer.write( " " );
        }
        writer.write( "</td>" );
        writer.write( "</tr>" );
        writer.write( ForesterUtil.LINE_SEPARATOR );
    }

    private static void writeValuesToTabWriter( final String species,
                                                final OntologizerResult ontologizer_result,
                                                final GoTerm got_term,
                                                final Writer writer ) throws IOException {
        writer.write( species );
        writer.write( "\t" );
        writer.write( got_term.getName() );
        writer.write( "\t" );
        writer.write( ontologizer_result.getGoId().getId() );
        writer.write( "\t" );
        writer.write( String.valueOf( ontologizer_result.getPAdjusted() ) );
        writer.write( "\t" );
        writer.write( String.valueOf( ontologizer_result.getP() ) );
        writer.write( "\t" );
        writer.write( String.valueOf( ontologizer_result.getPopTotal() ) );
        writer.write( "\t" );
        writer.write( String.valueOf( ontologizer_result.getPopTerm() ) );
        writer.write( "\t" );
        writer.write( String.valueOf( ontologizer_result.getStudyTotal() ) );
        writer.write( "\t" );
        writer.write( String.valueOf( ontologizer_result.getStudyTerm() ) );
        writer.write( "\t" );
        writer.write( String.valueOf( ontologizer_result.isTrivial() ) );
        writer.write( ForesterUtil.LINE_SEPARATOR );
    }
}
