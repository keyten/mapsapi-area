ymaps.modules.define('util.calculateArea', [], function (provide) {
    // Equatorial radius of Earth
    var RADIUS = 6378137;

    function calculateArea(feature) {
        var geoJsonGeometry = getGeoJsonGeometry(feature);
        return calculateJsonGeometryArea(geoJsonGeometry);
    }

    function getGeoJsonGeometry(feature) {
        if (feature.type === 'Feature') {
            return feature.geometry;
        } else if (feature.geometry && feature.geometry.getType) {
            if (feature.geometry.getType() === 'Circle') {
                return {
                    type: 'Circle',
                    coordinates: feature.geometry.getCoordinates(),
                    radius: feature.geometry.getRadius()
                };
            }
            return {
                type: feature.geometry.getType(),
                coordinates: feature.geometry.getCoordinates()
            };
        } else {
            throw new Error('util.calculateArea: Unknown input object.');
        }
    }

    function calculateJsonGeometryArea(geometry) {
        var area = 0;
        var i;
        switch (geometry.type) {
            case 'Polygon':
                if (isPolySelfIntersecting(geometry.coordinates)) {
                    return null;
                }

                return polygonArea(geometry.coordinates);
            case 'MultiPolygon':
                for (i = 0; i < geometry.coordinates.length; i++) {
                    area += polygonArea(geometry.coordinates[i]);
                }
                return area;
            case 'Rectangle':
                return polygonArea([[
                    geometry.coordinates[0],
                    [geometry.coordinates[0][0], geometry.coordinates[1][1]],
                    geometry.coordinates[1],
                    [geometry.coordinates[1][0], geometry.coordinates[0][1]],
                    geometry.coordinates[0]
                ]]);
            case 'Circle':
                return Math.PI * Math.pow(geometry.radius, 2);
            case 'Point':
            case 'MultiPoint':
            case 'LineString':
            case 'MultiLineString':
                return 0;
        }
    }

    function polygonArea(coords) {
        var area = 0;
        if (coords && coords.length > 0) {
            var signs = getRingsSigns(coords);
            for (var i = 0; i < coords.length; i++) {
                area += Math.abs(ringArea(coords[i])) * signs[i];
            }
        }
        return area;
    }

    /**
     * Modified version of https://github.com/mapbox/geojson-area
     * Calculate the approximate area of the polygon were it projected onto
     *     the earth.  Note that this area will be positive if ring is oriented
     *     clockwise, otherwise it will be negative.
     *
     * Reference:
     * Robert. G. Chamberlain and William H. Duquette, "Some Algorithms for
     *     Polygons on a Sphere", JPL Publication 07-03, Jet Propulsion
     *     Laboratory, Pasadena, CA, June 2007 https://trs.jpl.nasa.gov/handle/2014/40409
     *
     * Returns:
     * {Number} The approximate signed geodesic area of the polygon in square
     *     meters.
     */

    function ringArea(coords) {
        var p1;
        var p2;
        var p3;
        var lowerIndex;
        var middleIndex;
        var upperIndex;
        var i;
        var area = 0;
        var coordsLength = coords.length;
        var longitude = ymaps.meta.coordinatesOrder === 'latlong' ? 1 : 0;
        var latitude = ymaps.meta.coordinatesOrder === 'latlong' ? 0 : 1;

        if (coordsLength > 2) {
            for (i = 0; i < coordsLength; i++) {
                // i = N-2
                if (i === coordsLength - 2) {
                    lowerIndex = coordsLength - 2;
                    middleIndex = coordsLength - 1;
                    upperIndex = 0;
                } else if (i === coordsLength - 1) {
                    // i = N-1
                    lowerIndex = coordsLength - 1;
                    middleIndex = 0;
                    upperIndex = 1;
                } else {
                    // i = 0 to N-3
                    lowerIndex = i;
                    middleIndex = i + 1;
                    upperIndex = i + 2;
                }
                p1 = coords[lowerIndex];
                p2 = coords[middleIndex];
                p3 = coords[upperIndex];
                area += (rad(p3[longitude]) - rad(p1[longitude])) * Math.sin(rad(p2[latitude]));
            }

            area = area * RADIUS * RADIUS / 2;
        }

        return area;
    }

    function rad(_) {
        return _ * Math.PI / 180;
    }

    function isPolySelfIntersecting(rings) {
        for (var ri = 0; ri < rings.length; ri++) {
            for (var rj = ri; rj < rings.length; rj++) {
                var ring1 = rings[ri];
                var ring2 = rings[rj];
                for (var i = 0; i < ring1.length - 1; i++) {
                    for (var j = 0; j < ring2.length - 1; j++) {
                        if (ring1 === ring2 && (Math.abs(i - j) === 1 || Math.abs(i - j) === ring1.length - 2)) {
                            continue;
                        }

                        if (linesIntersect(ring1[i][0], ring1[i][1], ring1[i + 1][0], ring1[i + 1][1],
                            ring2[j][0], ring2[j][1], ring2[j + 1][0], ring2[j + 1][1])) {
                            return true;
                        }
                    }
                }
            }
        }
        return false;
    }

    function linesIntersect(ax0, ay0, ax1, ay1, bx0, by0, bx1, by1) {
        var d = (ax1 - ax0) * (by0 - by1) - (ay1 - ay0) * (bx0 - bx1);
        if (d === 0) {
            return false;
        }

        var t = (bx0 - ax0) * (by0 - by1) - (by0 - ay0) * (bx0 - bx1);
        var w = (ax1 - ax0) * (by0 - ay0) - (bx0 - ax0) * (ay1 - ay0);

        t /= d;
        w /= d;

        return t >= 0 && t <= 1 && w >= 0 && w <= 1;
    }

    function getRingsSigns(poly) {
        poly = geoCoords.geoPolyAsPixel(poly);
        poly.forEach(function (ring) {
            // isPointInRing thinks ring is not closed
            ring.pop();
        });

        var result = [];
        for (var i = 0; i < poly.length; i++) {
            result.push(getRingSign(poly, i));
        }
        return result;
    }

    var minLatDelta = 0.00000000000001;
    function getRingSign(poly, ringIndex) {
        var ring = poly[ringIndex];
        // looking for the bottomest vertex of the ring
        var lowestPoint = ring[0];
        ring.forEach(function (point) {
            if (point[1] > lowestPoint[1]) {
                lowestPoint = point;
            }
        });

        // the point a biiit bottomer (exactly not inside our ring but not touching another ring)
        lowestPoint[1] += minLatDelta;

        // считаем контуры, внутри которых оказалась точка
        var result = 1;
        poly.forEach(function (ring) {
            if (isPointInRing(lowestPoint[0], lowestPoint[1], ring)) {
                result = -result;
            }
        });
        return result;
    }

    function isPointInRing(x, y, ring) {
        var i = 0;
        var j = ring.length - 1;
        var c = false;
        for (; i < ring.length; i++) {
            if ((ring[i][1] > y) !== (ring[j][1] > y) &&
                (x < (ring[j][0] - ring[i][0]) * (y - ring[i][1]) / (ring[j][1] - ring[i][1]) + ring[i][0])) {
                c = !c;
            }
            j = i;
        }
        return c;
    }

    provide(calculateArea);
});
