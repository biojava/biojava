
package org.forester.surfacing;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import org.forester.go.GoId;
import org.forester.surfacing.BinaryDomainCombination.DomainCombinationType;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterUtil;

public class BasicGenomeWideCombinableDomains implements GenomeWideCombinableDomains {

    private final static NumberFormat                    FORMATTER                                  = new DecimalFormat( "0.0E0" );
    private static final Comparator<CombinableDomains>   DESCENDING_KEY_DOMAIN_COUNT_ORDER          = new Comparator<CombinableDomains>() {

                                                                                                        public int compare( final CombinableDomains d1,
                                                                                                                            final CombinableDomains d2 ) {
                                                                                                            if ( d1
                                                                                                                    .getKeyDomainCount() < d2
                                                                                                                    .getKeyDomainCount() ) {
                                                                                                                return 1;
                                                                                                            }
                                                                                                            else if ( d1
                                                                                                                    .getKeyDomainCount() > d2
                                                                                                                    .getKeyDomainCount() ) {
                                                                                                                return -1;
                                                                                                            }
                                                                                                            else {
                                                                                                                return d1
                                                                                                                        .getKeyDomain()
                                                                                                                        .getId()
                                                                                                                        .compareTo( d2
                                                                                                                                .getKeyDomain()
                                                                                                                                .getId() );
                                                                                                            }
                                                                                                        }
                                                                                                    };
    private static final Comparator<CombinableDomains>   DESCENDING_KEY_DOMAIN_PROTEINS_COUNT_ORDER = new Comparator<CombinableDomains>() {

                                                                                                        public int compare( final CombinableDomains d1,
                                                                                                                            final CombinableDomains d2 ) {
                                                                                                            if ( d1
                                                                                                                    .getKeyDomainProteinsCount() < d2
                                                                                                                    .getKeyDomainProteinsCount() ) {
                                                                                                                return 1;
                                                                                                            }
                                                                                                            else if ( d1
                                                                                                                    .getKeyDomainProteinsCount() > d2
                                                                                                                    .getKeyDomainProteinsCount() ) {
                                                                                                                return -1;
                                                                                                            }
                                                                                                            else {
                                                                                                                return d1
                                                                                                                        .getKeyDomain()
                                                                                                                        .getId()
                                                                                                                        .compareTo( d2
                                                                                                                                .getKeyDomain()
                                                                                                                                .getId() );
                                                                                                            }
                                                                                                        }
                                                                                                    };
    private static final Comparator<CombinableDomains>   DESCENDING_COMBINATIONS_COUNT_ORDER        = new Comparator<CombinableDomains>() {

                                                                                                        public int compare( final CombinableDomains d1,
                                                                                                                            final CombinableDomains d2 ) {
                                                                                                            if ( d1
                                                                                                                    .getNumberOfCombinableDomains() < d2
                                                                                                                    .getNumberOfCombinableDomains() ) {
                                                                                                                return 1;
                                                                                                            }
                                                                                                            else if ( d1
                                                                                                                    .getNumberOfCombinableDomains() > d2
                                                                                                                    .getNumberOfCombinableDomains() ) {
                                                                                                                return -1;
                                                                                                            }
                                                                                                            else {
                                                                                                                return d1
                                                                                                                        .getKeyDomain()
                                                                                                                        .getId()
                                                                                                                        .compareTo( d2
                                                                                                                                .getKeyDomain()
                                                                                                                                .getId() );
                                                                                                            }
                                                                                                        }
                                                                                                    };
    final private SortedMap<DomainId, CombinableDomains> _combinable_domains_map;
    final private Species                                _species;
    final private DomainCombinationType                  _dc_type;

    private BasicGenomeWideCombinableDomains( final Species species, final DomainCombinationType dc_type ) {
        _combinable_domains_map = new TreeMap<DomainId, CombinableDomains>();
        _species = species;
        _dc_type = dc_type;
    }

    private void add( final DomainId key, final CombinableDomains cdc ) {
        _combinable_domains_map.put( key, cdc );
    }

    public boolean contains( final DomainId key_id ) {
        return _combinable_domains_map.containsKey( key_id );
    }

    public CombinableDomains get( final DomainId key_id ) {
        return _combinable_domains_map.get( key_id );
    }

    public SortedMap<DomainId, CombinableDomains> getAllCombinableDomainsIds() {
        return _combinable_domains_map;
    }

    @Override
    public SortedSet<DomainId> getAllDomainIds() {
        final SortedSet<DomainId> domains = new TreeSet<DomainId>();
        for( final DomainId key : getAllCombinableDomainsIds().keySet() ) {
            final CombinableDomains cb = getAllCombinableDomainsIds().get( key );
            final List<DomainId> ds = cb.getAllDomains();
            for( final DomainId d : ds ) {
                domains.add( d );
            }
        }
        return domains;
    }

    @Override
    public DomainCombinationType getDomainCombinationType() {
        return _dc_type;
    }

    @Override
    public SortedSet<DomainId> getMostPromiscuosDomain() {
        final SortedSet<DomainId> doms = new TreeSet<DomainId>();
        final int max = ( int ) getPerGenomeDomainPromiscuityStatistics().getMax();
        for( final DomainId key : getAllCombinableDomainsIds().keySet() ) {
            final CombinableDomains cb = getAllCombinableDomainsIds().get( key );
            if ( cb.getNumberOfCombinableDomains() == max ) {
                doms.add( key );
            }
        }
        return doms;
    }

    @Override
    public DescriptiveStatistics getPerGenomeDomainPromiscuityStatistics() {
        final DescriptiveStatistics stats = new BasicDescriptiveStatistics();
        for( final DomainId key : getAllCombinableDomainsIds().keySet() ) {
            final CombinableDomains cb = getAllCombinableDomainsIds().get( key );
            stats.addValue( cb.getNumberOfCombinableDomains() );
        }
        return stats;
    }

    public int getSize() {
        return _combinable_domains_map.size();
    }

    public Species getSpecies() {
        return _species;
    }

    @Override
    public SortedSet<BinaryDomainCombination> toBinaryDomainCombinations() {
        final SortedSet<BinaryDomainCombination> binary_combinations = new TreeSet<BinaryDomainCombination>();
        for( final DomainId key : getAllCombinableDomainsIds().keySet() ) {
            final CombinableDomains cb = getAllCombinableDomainsIds().get( key );
            for( final BinaryDomainCombination b : cb.toBinaryDomainCombinations() ) {
                binary_combinations.add( b );
            }
        }
        return binary_combinations;
    }

    @Override
    public String toString() {
        return toStringBuilder( GenomeWideCombinableDomainsSortOrder.ALPHABETICAL_KEY_ID ).toString();
    }

    // Produces something like: 
    // 2-oxoacid_dh      5       5       2       4.8E-67   Biotin_lipoyl [4], E3_binding [3]
    public StringBuilder toStringBuilder( final GenomeWideCombinableDomainsSortOrder sort_order ) {
        final StringBuilder sb = new StringBuilder();
        final List<CombinableDomains> combinable_domains = new ArrayList<CombinableDomains>();
        for( final DomainId key : getAllCombinableDomainsIds().keySet() ) {
            final CombinableDomains cb = getAllCombinableDomainsIds().get( key );
            combinable_domains.add( cb );
        }
        if ( sort_order == GenomeWideCombinableDomainsSortOrder.KEY_DOMAIN_COUNT ) {
            Collections.sort( combinable_domains, BasicGenomeWideCombinableDomains.DESCENDING_KEY_DOMAIN_COUNT_ORDER );
        }
        else if ( sort_order == GenomeWideCombinableDomainsSortOrder.KEY_DOMAIN_PROTEINS_COUNT ) {
            Collections.sort( combinable_domains,
                              BasicGenomeWideCombinableDomains.DESCENDING_KEY_DOMAIN_PROTEINS_COUNT_ORDER );
        }
        else if ( sort_order == GenomeWideCombinableDomainsSortOrder.COMBINATIONS_COUNT ) {
            Collections.sort( combinable_domains, BasicGenomeWideCombinableDomains.DESCENDING_COMBINATIONS_COUNT_ORDER );
        }
        for( final CombinableDomains cb : combinable_domains ) {
            sb.append( ForesterUtil.pad( new StringBuffer( cb.getKeyDomain().toString() ), 18, ' ', false ) );
            sb.append( ForesterUtil.pad( new StringBuffer( "" + cb.getKeyDomainCount() ), 8, ' ', false ) );
            sb.append( ForesterUtil.pad( new StringBuffer( "" + cb.getKeyDomainProteinsCount() ), 8, ' ', false ) );
            sb.append( ForesterUtil.pad( new StringBuffer( "" + cb.getNumberOfCombinableDomains() ), 8, ' ', false ) );
            sb
                    .append( ForesterUtil
                            .pad( new StringBuffer( ""
                                          + FORMATTER
                                                  .format( cb.getKeyDomainConfidenceDescriptiveStatistics().median() ) ),
                                  10,
                                  ' ',
                                  false ) );
            sb.append( cb.getCombiningDomainIdsAsStringBuilder() );
            sb.append( ForesterUtil.getLineSeparator() );
        }
        return sb;
    }

    private static void countDomains( final Map<DomainId, Integer> domain_counts,
                                      final Map<DomainId, Integer> domain_protein_counts,
                                      final Map<DomainId, DescriptiveStatistics> stats,
                                      final Set<DomainId> saw_c,
                                      final DomainId id_i,
                                      final double support ) {
        if ( domain_counts.containsKey( id_i ) ) {
            domain_counts.put( id_i, 1 + domain_counts.get( ( id_i ) ) );
            if ( !saw_c.contains( id_i ) ) {
                domain_protein_counts.put( id_i, 1 + domain_protein_counts.get( ( id_i ) ) );
            }
        }
        else {
            stats.put( id_i, new BasicDescriptiveStatistics() );
            domain_counts.put( id_i, 1 );
            domain_protein_counts.put( id_i, 1 );
        }
        stats.get( id_i ).addValue( support );
        saw_c.add( id_i );
    }

    public static BasicGenomeWideCombinableDomains createInstance( final List<Protein> protein_list,
                                                                   final boolean ignore_combination_with_same_domain,
                                                                   final Species species ) {
        return createInstance( protein_list,
                               ignore_combination_with_same_domain,
                               species,
                               null,
                               DomainCombinationType.BASIC );
    }

    public static BasicGenomeWideCombinableDomains createInstance( final List<Protein> protein_list,
                                                                   final boolean ignore_combination_with_same_domain,
                                                                   final Species species,
                                                                   final DomainCombinationType dc_type ) {
        return createInstance( protein_list, ignore_combination_with_same_domain, species, null, dc_type );
    }

    public static BasicGenomeWideCombinableDomains createInstance( final List<Protein> protein_list,
                                                                   final boolean ignore_combination_with_same_domain,
                                                                   final Species species,
                                                                   final Map<DomainId, List<GoId>> domain_id_to_go_ids_map,
                                                                   final DomainCombinationType dc_type ) {
        final BasicGenomeWideCombinableDomains instance = new BasicGenomeWideCombinableDomains( species, dc_type );
        final Map<DomainId, Integer> domain_counts = new HashMap<DomainId, Integer>();
        final Map<DomainId, Integer> domain_protein_counts = new HashMap<DomainId, Integer>();
        final Map<DomainId, DescriptiveStatistics> stats = new HashMap<DomainId, DescriptiveStatistics>();
        for( final Protein protein : protein_list ) {
            if ( !protein.getSpecies().equals( species ) ) {
                throw new IllegalArgumentException( "species (" + protein.getSpecies()
                        + ") does not match species of combinable domains collection (" + species + ")" );
            }
            final Set<DomainId> saw_i = new HashSet<DomainId>();
            final Set<DomainId> saw_c = new HashSet<DomainId>();
            for( int i = 0; i < protein.getProteinDomains().size(); ++i ) {
                final Domain pd_i = protein.getProteinDomain( i );
                final DomainId id_i = pd_i.getDomainId();
                final int current_start = pd_i.getFrom();
                BasicGenomeWideCombinableDomains.countDomains( domain_counts,
                                                               domain_protein_counts,
                                                               stats,
                                                               saw_c,
                                                               id_i,
                                                               pd_i.getPerSequenceEvalue() );
                if ( !saw_i.contains( id_i ) ) {
                    if ( dc_type == DomainCombinationType.BASIC ) {
                        saw_i.add( id_i );
                    }
                    CombinableDomains domain_combination = null;
                    if ( instance.contains( id_i ) ) {
                        domain_combination = instance.get( id_i );
                    }
                    else {
                        if ( dc_type == DomainCombinationType.DIRECTED_ADJACTANT ) {
                            domain_combination = new AdjactantDirectedCombinableDomains( pd_i.getDomainId(), species );
                        }
                        else if ( dc_type == DomainCombinationType.DIRECTED ) {
                            domain_combination = new DirectedCombinableDomains( pd_i.getDomainId(), species );
                        }
                        else {
                            domain_combination = new BasicCombinableDomains( pd_i.getDomainId(), species );
                        }
                        if ( ( domain_id_to_go_ids_map != null )
                                && domain_id_to_go_ids_map.containsKey( pd_i.getDomainId() ) ) {
                            final List<GoId> go_ids = domain_id_to_go_ids_map.get( pd_i.getDomainId() );
                            for( final GoId go_id : go_ids ) {
                                domain_combination.getKeyDomain().addGoId( go_id );
                            }
                        }
                        instance.add( id_i, domain_combination );
                    }
                    final Set<DomainId> saw_j = new HashSet<DomainId>();
                    if ( ignore_combination_with_same_domain ) {
                        saw_j.add( id_i );
                    }
                    Domain closest = null;
                    for( int j = 0; j < protein.getNumberOfProteinDomains(); ++j ) {
                        if ( ( dc_type != DomainCombinationType.BASIC )
                                && ( current_start >= protein.getProteinDomain( j ).getFrom() ) ) {
                            continue;
                        }
                        if ( i != j ) {
                            final DomainId id = protein.getProteinDomain( j ).getDomainId();
                            if ( !saw_j.contains( id ) ) {
                                saw_j.add( id );
                                if ( dc_type != DomainCombinationType.DIRECTED_ADJACTANT ) {
                                    domain_combination
                                            .addCombinableDomain( protein.getProteinDomain( j ).getDomainId() );
                                }
                                else {
                                    if ( closest == null ) {
                                        closest = protein.getProteinDomain( j );
                                    }
                                    else {
                                        if ( protein.getProteinDomain( j ).getFrom() < closest.getFrom() ) {
                                            closest = protein.getProteinDomain( j );
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if ( ( dc_type == DomainCombinationType.DIRECTED_ADJACTANT ) && ( closest != null ) ) {
                        domain_combination.addCombinableDomain( closest.getDomainId() );
                    }
                }
            }
        }
        for( final DomainId key_id : domain_counts.keySet() ) {
            instance.get( key_id ).setKeyDomainCount( domain_counts.get( key_id ) );
            instance.get( key_id ).setKeyDomainProteinsCount( domain_protein_counts.get( key_id ) );
            instance.get( key_id ).setKeyDomainConfidenceDescriptiveStatistics( stats.get( key_id ) );
        }
        return instance;
    }
}
