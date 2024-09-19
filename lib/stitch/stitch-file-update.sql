drop index field_id_blocks;
drop index blocks_timestamp;
drop index blocks_x_min;
drop index blocks_y_min;
drop index blocks_z_min;
drop index blocks_x_max;
drop index blocks_y_max;
drop index blocks_z_max;

create index blocks_dims_and_field on blocks (x_min, y_min, z_min, x_max, y_max, z_max, field_id);

