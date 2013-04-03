// $Id: BasicGoTerm.java,v 1.26 2009/11/17 03:51:34 cmzmasek Exp $
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

package org.forester.go;

import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;

import org.forester.phylogeny.data.PhylogenyData;
import org.forester.util.ForesterUtil;

public class BasicGoTerm implements GoTerm {

    private final GoId           _id;
    private final String         _name;
    private final boolean        _is_obsolete;
    private final GoNameSpace    _namespace;
    private String               _definition;
    private List<GoId>           _alt_ids;
    private List<GoId>           _super_go_ids;
    private List<GoXRef>         _go_xrefs;
    private List<GoSubset>       _go_subsets;
    private String               _comment;
    private List<GoRelationship> _go_relationships;

    public BasicGoTerm( final GoId id, final String name, final GoNameSpace namespace, final boolean is_obsolete ) {
        if ( ( id == null ) || ForesterUtil.isEmpty( name ) || ( namespace == null ) ) {
            throw new IllegalArgumentException( "attempt to create GO term with empty id, name, or namespace" );
        }
        _id = id;
        _name = name;
        _namespace = namespace;
        _is_obsolete = is_obsolete;
        init();
    }

    public BasicGoTerm( final String id, final String name, final String namespace, final boolean is_obsolete ) {
        if ( ForesterUtil.isEmpty( id ) || ForesterUtil.isEmpty( name ) || ForesterUtil.isEmpty( namespace ) ) {
            throw new IllegalArgumentException( "attempt to create GO term with empty id, name, or namespace" );
        }
        _id = new GoId( id );
        _name = name;
        _namespace = new GoNameSpace( namespace );
        _is_obsolete = is_obsolete;
        init();
    }

    public StringBuffer asSimpleText() {
        return new StringBuffer( getGoId().toString() );
    }

    public StringBuffer asText() {
        return new StringBuffer( toString() );
    }

    /**
     * Compares based on GO id.
     * 
     */
    public int compareTo( final GoTerm go_term ) {
        return getGoId().compareTo( go_term.getGoId() );
    }

    /**
     * Makes a shallow copy.
     * 
     * 
     */
    public PhylogenyData copy() {
        final BasicGoTerm gt = new BasicGoTerm( getGoId(), getName(), getGoNameSpace(), isObsolete() );
        gt.setGoXrefs( getGoXRefs() );
        gt.setGoSubsets( getGoSubsets() );
        gt.setSuperTerms( getSuperGoIds() );
        gt.setAltIds( getAltIds() );
        gt.setDefinition( getDefinition() );
        return gt;
    }

    /**
     * Return true if both GO id and namespace are equal.
     * 
     */
    @Override
    public boolean equals( final Object o ) {
        if ( this == o ) {
            return true;
        }
        else if ( o == null ) {
            throw new IllegalArgumentException( "attempt to check go term equality to null" );
        }
        else if ( o.getClass() != this.getClass() ) {
            throw new IllegalArgumentException( "attempt to check go term equality to " + o + " [" + o.getClass() + "]" );
        }
        else {
            final GoTerm gt = ( GoTerm ) o;
            return getGoNameSpace().equals( gt.getGoNameSpace() ) && getGoId().equals( gt.getGoId() );
        }
    }

    public List<GoId> getAltIds() {
        return _alt_ids;
    }

    @Override
    public String getComment() {
        return _comment;
    }

    @Override
    public String getDefinition() {
        return _definition;
    }

    public GoId getGoId() {
        return _id;
    }

    public GoNameSpace getGoNameSpace() {
        return _namespace;
    }

    @Override
    public List<GoRelationship> getGoRelationships() {
        return _go_relationships;
    }

    @Override
    public List<GoSubset> getGoSubsets() {
        return _go_subsets;
    }

    public List<GoXRef> getGoXRefs() {
        return _go_xrefs;
    }

    public String getName() {
        return _name;
    }

    public List<GoId> getSuperGoIds() {
        return _super_go_ids;
    }

    /**
     * Hashcode is based on hashcode of GO id.
     * 
     * 
     */
    @Override
    public int hashCode() {
        return getGoId().hashCode();
    }

    private void init() {
        setGoXrefs( new ArrayList<GoXRef>() );
        setSuperTerms( new ArrayList<GoId>() );
        setAltIds( new ArrayList<GoId>() );
        setGoRelationships( new ArrayList<GoRelationship>() );
        setGoSubsets( new ArrayList<GoSubset>() );
        setDefinition( "" );
        setComment( "" );
    }

    public boolean isEqual( final PhylogenyData go_term ) {
        return equals( go_term );
    }

    public boolean isObsolete() {
        return _is_obsolete;
    }

    private void setAltIds( final List<GoId> alt_ids ) {
        _alt_ids = alt_ids;
    }

    public void setComment( final String comment ) {
        _comment = comment;
    }

    public void setDefinition( final String definition ) {
        _definition = definition;
    }

    private void setGoRelationships( final List<GoRelationship> go_relationships ) {
        _go_relationships = go_relationships;
    }

    public void setGoSubsets( final List<GoSubset> go_subsets ) {
        _go_subsets = go_subsets;
    }

    private void setGoXrefs( final List<GoXRef> xrefs ) {
        _go_xrefs = xrefs;
    }

    private void setSuperTerms( final List<GoId> super_terms ) {
        _super_go_ids = super_terms;
    }

    public StringBuffer toNHX() {
        throw new UnsupportedOperationException();
    }

    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        throw new UnsupportedOperationException();
    }

    @Override
    public String toString() {
        final StringBuffer sb = new StringBuffer();
        sb.append( getGoId() );
        sb.append( ": " );
        sb.append( getName() );
        sb.append( " [" );
        sb.append( getGoNameSpace() );
        sb.append( "]" );
        if ( isObsolete() ) {
            sb.append( " [is obsolete]" );
        }
        return sb.toString();
    }
}
