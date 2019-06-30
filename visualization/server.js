var connect = require('connect');
var serveStatic = require('serve-static');
var parser = require('biojs-io-newick');
connect().use(serveStatic(__dirname)).listen(8080, function(){
    console.log('Server running on 8080...');
});