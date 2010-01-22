// $Id: HmmscanPerDomainTableParser.java,v 1.1 2009/11/06 03:04:55 cmzmasek Exp
// $
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

package org.forester.io.parsers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import org.forester.surfacing.BasicDomain;
import org.forester.surfacing.BasicProtein;
import org.forester.surfacing.Domain;
import org.forester.surfacing.DomainId;
import org.forester.surfacing.Protein;
import org.forester.surfacing.SurfacingUtil;
import org.forester.util.ForesterUtil;

public final class HmmscanPerDomainTableParser {

    private static final String           RETRO                       = "RETRO";
    private static final String           PHAGE                       = "PHAGE";
    private static final String           VIR                         = "VIR";
    private static final String           TRANSPOS                    = "TRANSPOS";
    private static final String           RV                          = "RV";
    private static final String           GAG                         = "GAG_";
    private static final String           HCV                         = "HCV_";
    private static final String           HERPES                      = "HERPES_";
    private static final String           BACULO                      = "BACULO_";
    private static final int              E_VALUE_MAXIMUM_DEFAULT     = -1;
    private static final ReturnType       RETURN_TYPE_DEFAULT         = ReturnType.UNORDERED_PROTEIN_DOMAIN_COLLECTION_PER_PROTEIN;
    private static final boolean          IGNORE_DUFS_DEFAULT         = false;
    private static final int              MAX_ALLOWED_OVERLAP_DEFAULT = -1;
    private final Set<DomainId>           _filter;
    private final FilterType              _filter_type;
    private final File                    _input_file;
    private final String                  _species;
    private double                        _e_value_maximum;
    private Map<String, Double>           _individual_score_cutoffs;
    private boolean                       _ignore_dufs;
    private boolean                       _ignore_virus_like_ids;
    private int                           _max_allowed_overlap;
    private boolean                       _ignore_engulfed_domains;
    private ReturnType                    _return_type;
    private int                           _proteins_encountered;
    private int                           _proteins_ignored_due_to_filter;
    private int                           _proteins_stored;
    private int                           _domains_encountered;
    private int                           _domains_ignored_due_to_duf;
    private int                           _domains_ignored_due_to_overlap;
    private int                           _domains_ignored_due_to_e_value;
    private int                           _domains_ignored_due_to_individual_score_cutoff;
    private int                           _domains_stored;
    private SortedSet<DomainId>           _domains_stored_set;
    private long                          _time;
    private int                           _domains_ignored_due_to_negative_domain_filter;
    private Map<String, Integer>          _domains_ignored_due_to_negative_domain_filter_counts_map;
    private int                           _domains_ignored_due_to_virus_like_id;
    private Map<String, Integer>          _domains_ignored_due_to_virus_like_id_counts_map;
    private final INDIVIDUAL_SCORE_CUTOFF _ind_cutoff;

    public HmmscanPerDomainTableParser( final File input_file,
                                        final String species,
                                        final INDIVIDUAL_SCORE_CUTOFF individual_cutoff_applies_to ) {
        _input_file = input_file;
        _species = species;
        _filter = null;
        _filter_type = FilterType.NONE;
        _ind_cutoff = individual_cutoff_applies_to;
        init();
    }

    public HmmscanPerDomainTableParser( final File input_file,
                                        final String species,
                                        final Set<DomainId> filter,
                                        final FilterType filter_type,
                                        final INDIVIDUAL_SCORE_CUTOFF individual_cutoff_applies_to ) {
        _input_file = input_file;
        _species = species;
        _filter = filter;
        _filter_type = filter_type;
        _ind_cutoff = individual_cutoff_applies_to;
        init();
    }

    private void actuallyAddProtein( final List<Protein> proteins, final Protein current_protein ) {
        final List<Domain> l = current_protein.getProteinDomains();
        for( final Domain d : l ) {
            getDomainsStoredSet().add( d.getDomainId() );
        }
        proteins.add( current_protein );
        ++_proteins_stored;
    }

    private void addProtein( final List<Protein> proteins, Protein current_protein ) {
        if ( ( getMaxAllowedOverlap() != HmmscanPerDomainTableParser.MAX_ALLOWED_OVERLAP_DEFAULT )
                || isIgnoreEngulfedDomains() ) {
            final int domains_count = current_protein.getNumberOfProteinDomains();
            current_protein = SurfacingUtil.removeOverlappingDomains( getMaxAllowedOverlap(),
                                                                      isIgnoreEngulfedDomains(),
                                                                      current_protein );
            final int domains_removed = domains_count - current_protein.getNumberOfProteinDomains();
            _domains_stored -= domains_removed;
            _domains_ignored_due_to_overlap += domains_removed;
        }
        if ( ( getFilterType() == FilterType.POSITIVE_PROTEIN ) || ( getFilterType() == FilterType.NEGATIVE_PROTEIN ) ) {
            final Set<DomainId> domain_ids_in_protein = new HashSet<DomainId>();
            for( final Domain d : current_protein.getProteinDomains() ) {
                domain_ids_in_protein.add( d.getDomainId() );
            }
            domain_ids_in_protein.retainAll( getFilter() );
            if ( getFilterType() == FilterType.POSITIVE_PROTEIN ) {
                if ( domain_ids_in_protein.size() > 0 ) {
                    actuallyAddProtein( proteins, current_protein );
                }
                else {
                    ++_proteins_ignored_due_to_filter;
                }
            }
            else {
                if ( domain_ids_in_protein.size() < 1 ) {
                    actuallyAddProtein( proteins, current_protein );
                }
                else {
                    ++_proteins_ignored_due_to_filter;
                }
            }
        }
        else {
            actuallyAddProtein( proteins, current_protein );
        }
    }

    public int getDomainsEncountered() {
        return _domains_encountered;
    }

    public int getDomainsIgnoredDueToDuf() {
        return _domains_ignored_due_to_duf;
    }

    public int getDomainsIgnoredDueToEval() {
        return _domains_ignored_due_to_e_value;
    }

    public int getDomainsIgnoredDueToIndividualScoreCutoff() {
        return _domains_ignored_due_to_individual_score_cutoff;
    }

    public int getDomainsIgnoredDueToNegativeDomainFilter() {
        return _domains_ignored_due_to_negative_domain_filter;
    }

    public Map<String, Integer> getDomainsIgnoredDueToNegativeDomainFilterCountsMap() {
        return _domains_ignored_due_to_negative_domain_filter_counts_map;
    }

    public int getDomainsIgnoredDueToOverlap() {
        return _domains_ignored_due_to_overlap;
    }

    public Map<String, Integer> getDomainsIgnoredDueToVirusLikeIdCountsMap() {
        return _domains_ignored_due_to_virus_like_id_counts_map;
    }

    public int getDomainsIgnoredDueToVirusLikeIds() {
        return _domains_ignored_due_to_virus_like_id;
    }

    public int getDomainsStored() {
        return _domains_stored;
    }

    public SortedSet<DomainId> getDomainsStoredSet() {
        return _domains_stored_set;
    }

    private double getEValueMaximum() {
        return _e_value_maximum;
    }

    private Set<DomainId> getFilter() {
        return _filter;
    }

    private FilterType getFilterType() {
        return _filter_type;
    }

    public INDIVIDUAL_SCORE_CUTOFF getIndividualCutoffAppliesTo() {
        return _ind_cutoff;
    }

    private Map<String, Double> getIndividualScoreCutoffs() {
        return _individual_score_cutoffs;
    }

    private File getInputFile() {
        return _input_file;
    }

    private int getMaxAllowedOverlap() {
        return _max_allowed_overlap;
    }

    public int getProteinsEncountered() {
        return _proteins_encountered;
    }

    public int getProteinsIgnoredDueToFilter() {
        return _proteins_ignored_due_to_filter;
    }

    public int getProteinsStored() {
        return _proteins_stored;
    }

    private ReturnType getReturnType() {
        return _return_type;
    }

    private String getSpecies() {
        return _species;
    }

    public long getTime() {
        return _time;
    }

    private void init() {
        _e_value_maximum = HmmscanPerDomainTableParser.E_VALUE_MAXIMUM_DEFAULT;
        setIgnoreDufs( HmmscanPerDomainTableParser.IGNORE_DUFS_DEFAULT );
        setReturnType( HmmscanPerDomainTableParser.RETURN_TYPE_DEFAULT );
        _max_allowed_overlap = HmmscanPerDomainTableParser.MAX_ALLOWED_OVERLAP_DEFAULT;
        setIndividualScoreCutoffs( null );
        setIgnoreEngulfedDomains( false );
        setIgnoreVirusLikeIds( false );
        intitCounts();
    }

    private void intitCounts() {
        setDomainsStoredSet( new TreeSet<DomainId>() );
        setDomainsEncountered( 0 );
        setProteinsEncountered( 0 );
        setProteinsIgnoredDueToFilter( 0 );
        setDomainsIgnoredDueToNegativeFilter( 0 );
        setDomainsIgnoredDueToDuf( 0 );
        setDomainsIgnoredDueToEval( 0 );
        setDomainsIgnoredDueToIndividualScoreCutoff( 0 );
        setDomainsIgnoredDueToVirusLikeId( 0 );
        setDomainsIgnoredDueToOverlap( 0 );
        setDomainsStored( 0 );
        setProteinsStored( 0 );
        setTime( 0 );
        setDomainsIgnoredDueToVirusLikeIdCountsMap( new TreeMap<String, Integer>() );
        setDomainsIgnoredDueToNegativeDomainFilterCountsMap( new TreeMap<String, Integer>() );
    }

    private boolean isIgnoreDufs() {
        return _ignore_dufs;
    }

    private boolean isIgnoreEngulfedDomains() {
        return _ignore_engulfed_domains;
    }

    private boolean isIgnoreVirusLikeIds() {
        return _ignore_virus_like_ids;
    }

    public List<Protein> parse() throws IOException {
        if ( ( getIndividualCutoffAppliesTo() != INDIVIDUAL_SCORE_CUTOFF.NONE )
                && ( ( getIndividualScoreCutoffs() == null ) || ( getIndividualScoreCutoffs().size() < 1 ) ) ) {
            throw new IllegalStateException( "attempt to use individual cuttoffs with having set them" );
        }
        intitCounts();
        final Set<String> prev_queries = new HashSet<String>();
        final String error = ForesterUtil.isReadableFile( getInputFile() );
        if ( !ForesterUtil.isEmpty( error ) ) {
            throw new IOException( error );
        }
        final BufferedReader br = new BufferedReader( new FileReader( getInputFile() ) );
        String line;
        final List<Protein> proteins = new ArrayList<Protein>();
        Protein current_protein = null;
        int line_number = 0;
        final long start_time = new Date().getTime();
        String prev_query = "";
        int prev_qlen = -1;
        while ( ( line = br.readLine() ) != null ) {
            line_number++;
            if ( ForesterUtil.isEmpty( line ) || line.startsWith( "#" ) ) {
                continue;
            }
            // 0                    1           2    3                      4           5      6        7      8      9  10  11        12        13     14    15      16  17      18  19      20  21  22      
            // #                                                                              --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
            // # target name        accession   tlen query name             accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
            // #------------------- ---------- -----   -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
            // Ion_trans            PF00520.24   201 jgi|Nemve1|7|gw.28.1.1 -           1604  6.3e-169  557.4  95.3   1   4   1.5e-41     3e-38  130.8  11.1     3   171   140   307   139   346 0.81 Ion transport protein
            // Ion_trans            PF00520.24   201 jgi|Nemve1|7|gw.28.1.1 -           1604  6.3e-169  557.4  95.3   2   4   9.1e-45   1.8e-41  141.3  13.1     4   200   479   664   476   665 0.97 Ion transport protein
            // Ion_trans            PF00520.24   201 jgi|Nemve1|7|gw.28.1.1 -           1604  6.3e-169  557.4  95.3   3   4   5.2e-45     1e-41  142.1  14.0     1   201   900  1117   900  1117 0.96 Ion transport protein
            // Ion_trans            PF00520.24   201 jgi|Nemve1|7|gw.28.1.1 -           1604  6.3e-169  557.4  95.3   4   4   9.2e-51   1.8e-47  160.9  11.3     1   201  1217  1423  1217  1423 0.97 Ion transport protein
            // PKD_channel          PF08016.5    426 jgi|Nemve1|7|gw.28.1.1 -           1604   5.9e-19   67.4  70.5   1   8   0.00053       1.1    7.3   0.4   220   264   142   191   134   200 0.73 Polycystin cation channel
            final String tokens[] = line.split( "\\s+" );
            final String target_id = tokens[ 0 ];
            final String target_acc = tokens[ 1 ];
            final int tlen = parseInt( tokens[ 2 ], line_number, "tlen" );
            final String query = tokens[ 3 ];
            final String query_acc = tokens[ 4 ];
            final int qlen = parseInt( tokens[ 5 ], line_number, "qlen" );
            final double fs_e_value = parseDouble( tokens[ 6 ], line_number, "E-value" );
            final double fs_score = parseDouble( tokens[ 7 ], line_number, "score" );
            final int domain_number = parseInt( tokens[ 9 ], line_number, "count" );
            final int total_domains = parseInt( tokens[ 10 ], line_number, "total" );
            final double c_e_value = parseDouble( tokens[ 11 ], line_number, "c-Evalue" );
            final double i_e_value = parseDouble( tokens[ 12 ], line_number, "i-Evalue" );
            final double domain_score = parseDouble( tokens[ 13 ], line_number, "score" );
            final int hmm_from = parseInt( tokens[ 15 ], line_number, "hmm from" );
            final int hmm_to = parseInt( tokens[ 16 ], line_number, "hmm to" );
            final int ali_from = parseInt( tokens[ 17 ], line_number, "ali from" );
            final int ali_to = parseInt( tokens[ 18 ], line_number, "ali to" );
            final int env_from = parseInt( tokens[ 19 ], line_number, "env from" );
            final int env_to = parseInt( tokens[ 20 ], line_number, "env to" );
            ++_domains_encountered;
            if ( !query.equals( prev_query ) || ( qlen != prev_qlen ) ) {
                if ( query.equals( prev_query ) ) {
                    throw new IOException( "more than one protein named [" + query + "]" + " lengths: " + qlen + ", "
                            + prev_qlen );
                }
                if ( prev_queries.contains( query ) ) {
                    throw new IOException( "more than one protein named [" + query + "]" );
                }
                prev_query = query;
                prev_qlen = qlen;
                prev_queries.add( query );
                if ( ( current_protein != null ) && ( current_protein.getProteinDomains().size() > 0 ) ) {
                    addProtein( proteins, current_protein );
                }
                if ( getReturnType() == ReturnType.UNORDERED_PROTEIN_DOMAIN_COLLECTION_PER_PROTEIN ) {
                    current_protein = new BasicProtein( query, getSpecies() );
                }
                else {
                    throw new IllegalArgumentException( "unknown return type" );
                }
            }
            boolean failed_cutoff = false;
            if ( getIndividualCutoffAppliesTo() != INDIVIDUAL_SCORE_CUTOFF.NONE ) {
                if ( getIndividualScoreCutoffs().containsKey( target_id ) ) {
                    final double cutoff = getIndividualScoreCutoffs().get( target_id );
                    if ( getIndividualCutoffAppliesTo() != INDIVIDUAL_SCORE_CUTOFF.FULL_SEQUENCE ) {
                        if ( fs_score < cutoff ) {
                            failed_cutoff = true;
                        }
                    }
                    else if ( getIndividualCutoffAppliesTo() != INDIVIDUAL_SCORE_CUTOFF.DOMAIN ) {
                        if ( domain_score < cutoff ) {
                            failed_cutoff = true;
                        }
                    }
                }
                else {
                    throw new IOException( "could not find a score cutoff value for domain id \"" + target_id
                            + "\" [line " + line_number + "] in [" + getInputFile().getCanonicalPath() + "]" );
                }
            }
            final String uc_id = target_id.toUpperCase();
            if ( failed_cutoff ) {
                ++_domains_ignored_due_to_individual_score_cutoff;
            }
            else if ( ali_from == ali_to ) {
                //Ignore
            }
            else if ( ( getEValueMaximum() != HmmscanPerDomainTableParser.E_VALUE_MAXIMUM_DEFAULT )
                    && ( fs_e_value > getEValueMaximum() ) ) {
                ++_domains_ignored_due_to_e_value;
            }
            else if ( isIgnoreDufs() && uc_id.startsWith( "DUF" ) ) {
                ++_domains_ignored_due_to_duf;
            }
            else if ( isIgnoreVirusLikeIds()
                    && ( uc_id.contains( VIR ) || uc_id.contains( PHAGE ) || uc_id.contains( RETRO )
                            || uc_id.contains( TRANSPOS ) || uc_id.startsWith( RV ) || uc_id.startsWith( GAG )
                            || uc_id.startsWith( HCV ) || uc_id.startsWith( HERPES ) || uc_id.startsWith( BACULO ) ) ) {
                ForesterUtil.increaseCountingMap( getDomainsIgnoredDueToVirusLikeIdCountsMap(), target_id );
                ++_domains_ignored_due_to_virus_like_id;
            }
            else if ( ( getFilterType() == FilterType.NEGATIVE_DOMAIN )
                    && getFilter().contains( new DomainId( target_id ) ) ) {
                ++_domains_ignored_due_to_negative_domain_filter;
                ForesterUtil.increaseCountingMap( getDomainsIgnoredDueToNegativeDomainFilterCountsMap(), target_id );
            }
            else {
                try {
                    final Domain pd = new BasicDomain( target_id,
                                                       ali_from,
                                                       ali_to,
                                                       ( short ) domain_number,
                                                       ( short ) total_domains,
                                                       fs_e_value,
                                                       fs_score,
                                                       i_e_value,
                                                       domain_score );
                    current_protein.addProteinDomain( pd );
                }
                catch ( final IllegalArgumentException e ) {
                    throw new IOException( "problem with domain parsing at line " + line_number + "[" + line + "]: "
                            + e.getMessage() );
                }
                ++_domains_stored;
            }
        } // while ( ( line = br.readLine() ) != null )
        if ( ( current_protein != null ) && ( current_protein.getProteinDomains().size() > 0 ) ) {
            addProtein( proteins, current_protein );
        }
        setProteinsEncountered( prev_queries.size() );
        setTime( new Date().getTime() - start_time );
        return proteins;
    }

    private double parseDouble( final String double_str, final int line_number, final String label ) throws IOException {
        double d = -1;
        try {
            d = Double.valueOf( double_str ).doubleValue();
        }
        catch ( final NumberFormatException e ) {
            throw new IOException( "could not parse \" +label + \" from \"" + double_str + "\" [line " + line_number
                    + "] in [" + getInputFile().getCanonicalPath() + "]" );
        }
        return d;
    }

    private int parseInt( final String double_str, final int line_number, final String label ) throws IOException {
        int i = -1;
        try {
            i = Integer.valueOf( double_str ).intValue();
        }
        catch ( final NumberFormatException e ) {
            throw new IOException( "could not parse \"" + label + "\" from \"" + double_str + "\" [line " + line_number
                    + "] in [" + getInputFile().getCanonicalPath() + "]" );
        }
        return i;
    }

    private void setDomainsEncountered( final int domains_encountered ) {
        _domains_encountered = domains_encountered;
    }

    private void setDomainsIgnoredDueToDuf( final int domains_ignored_due_to_duf ) {
        _domains_ignored_due_to_duf = domains_ignored_due_to_duf;
    }

    public void setDomainsIgnoredDueToEval( final int domains_ignored_due_to_e_value ) {
        _domains_ignored_due_to_e_value = domains_ignored_due_to_e_value;
    }

    public void setDomainsIgnoredDueToIndividualScoreCutoff( final int domains_ignored_due_to_individual_score_cutoff ) {
        _domains_ignored_due_to_individual_score_cutoff = domains_ignored_due_to_individual_score_cutoff;
    }

    private void setDomainsIgnoredDueToNegativeDomainFilterCountsMap( final Map<String, Integer> domains_ignored_due_to_negative_domain_filter_counts_map ) {
        _domains_ignored_due_to_negative_domain_filter_counts_map = domains_ignored_due_to_negative_domain_filter_counts_map;
    }

    private void setDomainsIgnoredDueToNegativeFilter( final int domains_ignored_due_to_negative_domain_filter ) {
        _domains_ignored_due_to_negative_domain_filter = domains_ignored_due_to_negative_domain_filter;
    }

    private void setDomainsIgnoredDueToOverlap( final int domains_ignored_due_to_overlap ) {
        _domains_ignored_due_to_overlap = domains_ignored_due_to_overlap;
    }

    private void setDomainsIgnoredDueToVirusLikeId( final int i ) {
        _domains_ignored_due_to_virus_like_id = i;
    }

    private void setDomainsIgnoredDueToVirusLikeIdCountsMap( final Map<String, Integer> domains_ignored_due_to_virus_like_id_counts_map ) {
        _domains_ignored_due_to_virus_like_id_counts_map = domains_ignored_due_to_virus_like_id_counts_map;
    }

    private void setDomainsStored( final int domains_stored ) {
        _domains_stored = domains_stored;
    }

    private void setDomainsStoredSet( final SortedSet<DomainId> _storeddomains_stored ) {
        _domains_stored_set = _storeddomains_stored;
    }

    public void setEValueMaximum( final double e_value_maximum ) {
        if ( e_value_maximum < 0.0 ) {
            throw new IllegalArgumentException( "attempt to set the maximum E-value to a negative value" );
        }
        _e_value_maximum = e_value_maximum;
    }

    public void setIgnoreDufs( final boolean ignore_dufs ) {
        _ignore_dufs = ignore_dufs;
    }

    /**
     * To ignore domains which are completely engulfed by domains (individual
     * ones or stretches of overlapping ones) with better support values.
     * 
     * 
     * @param ignored_engulfed_domains
     */
    public void setIgnoreEngulfedDomains( final boolean ignore_engulfed_domains ) {
        _ignore_engulfed_domains = ignore_engulfed_domains;
    }

    public void setIgnoreVirusLikeIds( final boolean ignore_virus_like_ids ) {
        _ignore_virus_like_ids = ignore_virus_like_ids;
    }

    /**
     * Sets the individual  score cutoff values (for example, gathering
     * thresholds from Pfam). Domain ids are the keys, cutoffs the values.
     * 
     * @param individual_score_cutoffs
     */
    public void setIndividualScoreCutoffs( final Map<String, Double> individual_score_cutoffs ) {
        _individual_score_cutoffs = individual_score_cutoffs;
    }

    public void setMaxAllowedOverlap( final int max_allowed_overlap ) {
        if ( max_allowed_overlap < 0 ) {
            throw new IllegalArgumentException( "Attempt to set max allowed overlap to less than zero." );
        }
        _max_allowed_overlap = max_allowed_overlap;
    }

    private void setProteinsEncountered( final int proteins_encountered ) {
        _proteins_encountered = proteins_encountered;
    }

    private void setProteinsIgnoredDueToFilter( final int proteins_ignored_due_to_filter ) {
        _proteins_ignored_due_to_filter = proteins_ignored_due_to_filter;
    }

    private void setProteinsStored( final int proteins_stored ) {
        _proteins_stored = proteins_stored;
    }

    public void setReturnType( final ReturnType return_type ) {
        _return_type = return_type;
    }

    private void setTime( final long time ) {
        _time = time;
    }

    public static enum FilterType {
        NONE, POSITIVE_PROTEIN, NEGATIVE_PROTEIN, NEGATIVE_DOMAIN
    }

    static public enum INDIVIDUAL_SCORE_CUTOFF {
        FULL_SEQUENCE, DOMAIN, NONE;
    }

    public static enum ReturnType {
        UNORDERED_PROTEIN_DOMAIN_COLLECTION_PER_PROTEIN
    }
}
