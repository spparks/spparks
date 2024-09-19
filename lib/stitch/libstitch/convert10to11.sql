--
-- convert a Stitch-IO 1.0 database to the new 1.1 format with the blobs separate from the block metadata
--
-- create a view to make the rowid visible for our foreign key
-- create the new blobs table using that view;
create view blocks_split as select rowid as id, state from blocks;
create table blocks_blobs as select id, state from blocks_split;

-- alter table blocks drop column state; has to be done as below due to SQLite limitations
PRAGMA foreign_keys=off;

BEGIN TRANSACTION;

create table new_blocks (field_id int not null, timestamp int not null, x_min int not null, y_min int not null, z_min int not null, x_max int not null, y_max int not null, z_max int not null, block_id int not null);

insert into new_blocks (field_id, timestamp, x_min, y_min, z_min, x_max, y_max, z_max, block_id)
select field_id, timestamp, x_min, y_min, z_min, x_max, y_max, z_max, rowid from blocks;

drop view blocks_split;
drop table blocks;

alter table new_blocks rename to blocks;

commit;

PRAGMA foreign_keys=on;

-- recreate the old index structures and add two new ones
create index field_id_blocks on blocks (field_id);
create index blocks_timestamp on blocks (timestamp);
create index blocks_x_min on blocks (x_min);
create index blocks_y_min on blocks (y_min);
create index blocks_z_min on blocks (z_min);
create index blocks_x_max on blocks (x_max);
create index blocks_y_max on blocks (y_max);
create index blocks_z_max on blocks (z_max);
create index blocks_block_id on blocks (block_id);
create index blocks_blobs_id on blocks_blobs (id);

--
-- change the globals table to inlcude a database version field (oversight for version 1.0)
--
alter table globals add stitch_version integer;
update globals set stitch_version=11;

-- clean out empty database space to reduce size;
vacuum;
