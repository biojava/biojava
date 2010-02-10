// $Id: PropertiesMap.java,v 1.4 2009/06/09 00:05:05 cmzmasek Exp $
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

package org.forester.phylogeny.data;

import java.io.IOException;
import java.io.Writer;
import java.util.SortedMap;
import java.util.TreeMap;

public class PropertiesMap implements PhylogenyData {

    private final SortedMap<String, Property> _properties;

    public PropertiesMap() {
        _properties = new TreeMap<String, Property>();
    }

    public void addProperty( final Property property ) throws IllegalArgumentException {
        if ( getProperties().containsKey( property.getRef() ) ) {
            throw new IllegalArgumentException( "ref [" + property.getRef() + "] is already present" );
        }
        getProperties().put( property.getRef(), property );
    }

    @Override
    public StringBuffer asSimpleText() {
        final StringBuffer sb = new StringBuffer();
        boolean first = true;
        for( final String ref : getPropertyRefs() ) {
            if ( first ) {
                first = false;
            }
            else {
                sb.append( " " );
            }
            sb.append( getProperty( ref ).asText() );
        }
        return sb;
    }

    @Override
    public StringBuffer asText() {
        return asSimpleText();
    }

    @Override
    public PhylogenyData copy() {
        final PropertiesMap new_one = new PropertiesMap();
        for( final String r : getProperties().keySet() ) {
            new_one.addProperty( getProperties().get( r ) );
        }
        return new_one;
    }

    public SortedMap<String, Property> getProperties() {
        return _properties;
    }

    public Property[] getPropertiesArray() {
        final Property[] a = new Property[ getProperties().size() ];
        int i = 0;
        for( final String r : getProperties().keySet() ) {
            a[ i++ ] = getProperties().get( r );
        }
        return a;
    }

    public Property getProperty( final String ref ) throws IllegalArgumentException {
        Property p = null;
        if ( getProperties() != null ) {
            p = getProperties().get( ref );
        }
        if ( p != null ) {
            return p;
        }
        else {
            throw new IllegalArgumentException( "Ref [" + ref + "] is not present" );
        }
    }

    /**
     * Returns all property refs of this PhylogenyNode as String array.
     */
    public String[] getPropertyRefs() {
        if ( getProperties() == null ) {
            return new String[ 0 ];
        }
        final Property[] properties = getPropertiesArray();
        final String[] refs = new String[ properties.length ];
        for( int i = 0; i < properties.length; ++i ) {
            refs[ i ] = properties[ i ].getRef();
        }
        return refs;
    }

    @Override
    public boolean isEqual( final PhylogenyData data ) {
        throw new UnsupportedOperationException();
    }

    public boolean refExists( final String ref ) {
        if ( getProperties() != null ) {
            for( final String r : getProperties().keySet() ) {
                if ( r.equalsIgnoreCase( ref ) ) {
                    return true;
                }
            }
        }
        return false;
    }

    @Override
    public StringBuffer toNHX() {
        final StringBuffer sb = new StringBuffer();
        if ( getProperties() != null ) {
            for( final String ref : getProperties().keySet() ) {
                sb.append( getProperties().get( ref ).toNHX() );
            }
        }
        return sb;
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        if ( getProperties() != null ) {
            for( final String ref : getProperties().keySet() ) {
                getProperties().get( ref ).toPhyloXML( writer, level, indentation );
            }
        }
    }

    @Override
    public String toString() {
        return asSimpleText().toString();
    }
}
